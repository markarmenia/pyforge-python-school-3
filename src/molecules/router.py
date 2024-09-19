from fastapi import APIRouter, HTTPException, Query
from typing import List, Optional
from celery.result import AsyncResult
from src.celery_worker import celery_app
from src.tasks import substructure_search_task
from rdkit import Chem
import logging
import json
from src.molecules.dao import MoleculeDAO
from src.molecules.schema import MoleculeResponse, MoleculeAdd, MoleculeUpdate
from src.cache import redis_client
from hashlib import sha256
from contextlib import asynccontextmanager

logger = logging.getLogger(__name__)

router = APIRouter()

@asynccontextmanager
async def get_redis_client():
    client = await redis_client()
    try:
        yield client
    finally:
        await client.aclose()

def get_cache_key(*args, **kwargs) -> str:
    key = json.dumps({"args": args, "kwargs": kwargs}, sort_keys=True)
    return sha256(key.encode()).hexdigest()

async def invalidate_molecule_cache(redis_client_instance):
    # Delete all 'get_all*' keys
    keys_to_delete_all = await redis_client_instance.keys("get_all*")
    if keys_to_delete_all:
        await redis_client_instance.delete(*keys_to_delete_all)
    
    # Delete all 'molecule_by_id:*' keys
    keys_to_delete_by_id = await redis_client_instance.keys("molecule_by_id:*")
    if keys_to_delete_by_id:
        await redis_client_instance.delete(*keys_to_delete_by_id)
    
    logger.info("Molecule cache invalidated for both list and individual retrievals")
# Add molecule (smiles) and its identifier
@router.post("/molecules", status_code=201, tags=["Molecules"], summary="Add new molecules to the DB", response_description="Molecule added successfully")
async def add_molecule(molecule: MoleculeAdd) -> dict:
    logger.info(f"Received request to add molecule with SMILES: {molecule.smiles}")

    substructure = Chem.MolFromSmiles(molecule.smiles)
    if not substructure:
        logger.warning(f"Invalid SMILES molecule received: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")

    molecule_data = molecule.model_dump()
    try:
        new_molecule_id = await MoleculeDAO.add_molecule(**molecule_data)
        logger.info(f"Molecule added successfully with ID: {new_molecule_id}")
        
        # Invalidate the cache
        async with get_redis_client() as redis_client_instance:
            await invalidate_molecule_cache(redis_client_instance)
        
        return {"message": "The molecule is added!", "molecule": molecule}
    except HTTPException as e:
        logger.error(f"HTTPException occurred: {e.detail}")
        raise e
    except Exception as e:
        logger.error(f"Unexpected error while adding molecule: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Error adding the molecule")

@router.get("/molecules/{molecule_id}", tags=["Molecules"], summary="Retrieve molecule by ID", response_description="Molecule retrieved successfully")
async def get_molecule_by_id(molecule_id: int):
    cache_key = f"molecule_by_id:{molecule_id}"
    
    async with get_redis_client() as redis_client_instance:
        cached_result = await redis_client_instance.get(cache_key)
        if cached_result:
            logger.info(f"Cache hit for key: {cache_key}")
            return json.loads(cached_result)

    logger.info(f"Received request to retrieve molecule with ID: {molecule_id}")
    rez = await MoleculeDAO.find_full_data(molecule_id=molecule_id)
    if rez is None:
        logger.error(f"Molecule with ID {molecule_id} not found")
        raise HTTPException(status_code=404, detail=f"Molecule with id {molecule_id} does not exist!")

    # Convert Pydantic model to dict for serialization
    rez_dict = rez.model_dump()
    logger.info(f"Molecule with ID {molecule_id} retrieved successfully")
    
    async with get_redis_client() as redis_client_instance:
        await redis_client_instance.setex(cache_key, 300, json.dumps(rez_dict))
    
    return rez_dict

# Updating a molecule by identifier
@router.put("/molecules/{molecule_id}", tags=["Molecules"], summary="Update molecule by ID", response_model=MoleculeResponse)
async def update_molecule(molecule_id: int, updated_molecule: MoleculeUpdate):
    logger.info(f"Received request to update molecule with ID {molecule_id} with new data: {updated_molecule.smiles}")

    substructure = Chem.MolFromSmiles(updated_molecule.smiles)
    if not substructure:
        logger.warning(f"Invalid SMILES molecule received for update: {updated_molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")

    try:
        updated = await MoleculeDAO.update(molecule_id, updated_molecule)
        if updated is None:
            logger.warning(f"Molecule with ID {molecule_id} not found for update")
            raise HTTPException(status_code=404, detail="Molecule not found")
        logger.info(f"Molecule with ID {molecule_id} updated successfully")
        
        # Invalidate the cache
        async with get_redis_client() as redis_client_instance:
            await invalidate_molecule_cache(redis_client_instance)
        
        return updated
    except HTTPException as e:
        logger.error(f"HTTPException occurred during update: {e.detail}")
        raise e
    except Exception as e:
        logger.error(f"Unexpected error while updating molecule with ID {molecule_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Unexpected error occurred")

# Delete a molecule by identifier
@router.delete("/molecules/{molecule_id}", tags=["Molecules"], summary="Delete molecule by ID", response_description="Molecule Deleted")
async def delete_molecule(molecule_id: int) -> dict:
    logger.info(f"Received request to delete molecule with ID {molecule_id}")
    check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
    if check:
        logger.info(f"Molecule with ID {molecule_id} deleted successfully")
        
        # Invalidate the cache
        async with get_redis_client() as redis_client_instance:
            await invalidate_molecule_cache(redis_client_instance)
        
        return {"message": f"The molecule with id {molecule_id} is deleted!"}
    else:
        logger.error(f"Molecule with ID {molecule_id} not found for deletion")
        raise HTTPException(status_code=404, detail=f"Molecule with id {molecule_id} does not exist!")

# List all molecules
@router.get("/molecules", tags=["Search"], summary="Retrieve all molecules", response_description="Molecules retrieved successfully")
async def get_all(limit: Optional[int] = Query(None, description="Limit the number of molecules to retrieve")) -> List[MoleculeResponse]:
    cache_key = f"get_all:{limit}" if limit else "get_all"
    
    async with get_redis_client() as redis_client_instance:
        cached_results = await redis_client_instance.get(cache_key)
        if cached_results:
            logger.info(f"Cache hit for key: {cache_key}")
            return json.loads(cached_results)

    if limit:
        logger.info(f"Received request to retrieve {limit} molecules")
    else:
        logger.info(f"Received request to retrieve ALL molecules")
    
    try:
        results = await MoleculeDAO.find_all_molecules(limit)
        logger.info(f"{len(results)} molecules retrieved successfully")
        
        # Cache the results
        async with get_redis_client() as redis_client_instance:
            await redis_client_instance.setex(cache_key, 300, json.dumps(results))
        
        return results
    except Exception as e:
        logger.error(f"Unexpected error during search: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Error performing search")

# Substructure search for all added molecules
@router.get("/substructures", tags=["Search"], summary="Start substructure search", response_description="Task started")
async def start_substructure_search(smiles: str):
    logger.info(f"Received request to start substructure search for SMILES: {smiles}")
    try:
        task = substructure_search_task.delay(smiles)  # Submit task to Celery
        return {"task_id": task.id, "message": "Task started, use /substructures/{task_id} to get the results"}
    except Exception as e:
        logger.error(f"Error starting substructure search task: {e}")
        raise HTTPException(status_code=500, detail="Error starting substructure search task")

# Get the substructure search results
@router.get("/substructures/{task_id}", tags=["Search"], summary="Get substructure search results", response_description="Search results")
async def get_substructure_search_results(task_id: str):
    logger.info(f"Received request for substructure search results for task_id: {task_id}")
    
    task = AsyncResult(task_id, app=celery_app)
    if task.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task.state == 'SUCCESS':
        result = task.result
        return {"task_id": task_id, "status": "Task completed", "result": result}
    elif task.state == 'FAILURE':
        return {"task_id": task_id, "status": "Task failed", "error": str(task.info)}
    else:
        return {"task_id": task_id, "status": "Unknown task state"}