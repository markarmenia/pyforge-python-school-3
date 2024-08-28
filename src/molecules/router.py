import logging
from typing import List, Optional
from fastapi import APIRouter, HTTPException, Query
from rdkit import Chem
from src.molecules.dao import MoleculeDAO
from src.molecules.schema import MoleculeResponse, MoleculeAdd, MoleculeUpdate

logger = logging.getLogger(__name__)

router = APIRouter()

# Add molecule (smiles) and its identifier
@router.post("/molecules", status_code=201, tags=["Molecules"], summary="Add new molecules to the DB", response_description="Molecule added successfully")
async def add_molecule(molecule: MoleculeAdd) -> dict:
    '''
    Create a molecule with the following details:

    **smiles**: Molecule in SMILES format  
    '''
    logger.info(f"Received request to add molecule with SMILES: {molecule.smiles}")

    substructure = Chem.MolFromSmiles(molecule.smiles)
    if not substructure:
        logger.warning(f"Invalid SMILES molecule received: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")

    molecule_data = molecule.model_dump()
    try:
        new_molecule_id = await MoleculeDAO.add_molecule(**molecule_data)
        logger.info(f"Molecule added successfully with ID: {new_molecule_id}")
        return {"message": "The molecule is added!", "molecule": molecule}
    except HTTPException as e:
        logger.error(f"HTTPException occurred: {e.detail}")
        raise e
    except Exception as e:
        logger.error(f"Unexpected error while adding molecule: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Error adding the molecule")

# Get molecule by identifier
@router.get("/molecules/{molecule_id}", tags=["Molecules"], summary="Retrieve molecule by ID", response_description="Molecule retrieved successfully")
async def get_molecule_by_id(molecule_id: int):
    '''
    Get molecule by identifier:

    **id**: Required
    '''
    logger.info(f"Received request to retrieve molecule with ID: {molecule_id}")
    try:
        rez = await MoleculeDAO.find_full_data(molecule_id=molecule_id)
        if rez is None:
            logger.warning(f"Molecule with ID {molecule_id} not found")
            raise HTTPException(status_code=404, detail=f"Molecule with id {molecule_id} does not exist!")
        logger.info(f"Molecule with ID {molecule_id} retrieved successfully")
        return rez
    except Exception as e:
        logger.error(f"Unexpected error while retrieving molecule with ID {molecule_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Error retrieving the molecule")

# Updating a molecule by identifier
@router.put("/molecules/{molecule_id}", tags=["Molecules"], summary="Update molecule by ID", response_model=MoleculeResponse)
async def update_molecule(molecule_id: int, updated_molecule: MoleculeUpdate):
    '''
    Updating a molecule by identifier:

    **id**: Required
    **smiles**: Molecule in SMILES format 
    '''
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
    '''
    Delete a molecule by identifier

    **id**: Required
    '''
    logger.info(f"Received request to delete molecule with ID {molecule_id}")
    try:
        check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
        if check:
            logger.info(f"Molecule with ID {molecule_id} deleted successfully")
            return {"message": f"The molecule with id {molecule_id} is deleted!"}
        else:
            logger.warning(f"Molecule with ID {molecule_id} not found for deletion")
            raise HTTPException(status_code=404, detail=f"Molecule with id {molecule_id} does not exist!")
    except Exception as e:
        logger.error(f"Unexpected error while deleting molecule with ID {molecule_id}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Error deleting the molecule")

# List all molecules
@router.get("/molecules", tags=["Search"], summary="Retrieve all molecules", response_description="Molecules retrieved successfully")
async def get_all(limit: Optional[int] = Query(None, description="Limit the number of molecules to retrieve")) -> List[MoleculeResponse]:
    '''
    List all molecules 
    **Optional**: Limit the required number of molecules in the response. 

    '''
    if limit:
        logger.info(f"Received request to retrieve {limit} molecules")
    else:
        logger.info(f"Received request to retrieve ALL molecules")
    try:
        results= [molecule async for molecule in MoleculeDAO.find_all_molecules(limit)]
        logger.info(f" {len(results)} molecules retrieved successfully")
        return results
    except Exception as e:
        logger.error(f"Unexpected error during search: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Error performing search")


# Substructure search for all added molecules
@router.get("/substructures", tags=["Search"], summary="Substructure search molecules", response_description="Substructure match molecules retrieved successfully")
async def substructure_search(smiles: str):
    '''
    Substructure search

    **substructure**: SMILES molecule required
    '''
    logger.info(f"Received request for substructure search with SMILES: {smiles}")
    try:
        matches = await MoleculeDAO.substructure_search(smiles)
        if not matches:
            logger.warning("No matching molecules found for substructure search")
            raise HTTPException(status_code=404, detail="No matching molecules found")
        logger.info(f"Substructure search found {len(matches)} matches")
        return matches
    except ValueError as e:
        logger.error(f"ValueError during substructure search: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Unexpected error during substructure search: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Error performing substructure search")
