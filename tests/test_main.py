import json
import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient
from unittest.mock import patch, AsyncMock, MagicMock
from src.molecules.router import router
from src.molecules.dao import MoleculeDAO
from src.molecules.schema import MoleculeResponse, MoleculeUpdate


client = TestClient(router)

#Test 1 adding molecule (smiles)
@patch.object(MoleculeDAO, 'add_molecule', new_callable=AsyncMock)
def test_add_molecule(mock_add_molecule):
    mock_add_molecule.return_value = 1
    response = client.post("/molecules", json={"smiles": "c1ccccc1"})
    
    assert response.status_code == 201
    assert response.json() == {"message": "The molecule is added!", "molecule": {"smiles": "c1ccccc1"}}
    mock_add_molecule.assert_called_once()

#Test 2 adding molecule with an invalid smiles structure
@patch.object(MoleculeDAO, 'add_molecule', new_callable=AsyncMock)
def test_add_invalid_molecule(mock_MolFromSmiles):
    mock_MolFromSmiles.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.post("/molecules", json={"smiles": "test"})
    
    assert excinfo.value.status_code == 400
    assert excinfo.value.detail == "Invalid SMILES molecule"

# Test 3 Adding molecule with a SMILES that already exists
@patch.object(MoleculeDAO, 'add_molecule', new_callable=AsyncMock)
def test_add_existing_molecule(mock_add_molecule):
    mock_add_molecule.side_effect = [None, HTTPException(status_code=409, detail="Molecule already exists")]
    response = client.post("/molecules", json={"smiles": "c1ccccc1"})
    assert response.status_code == 201
    with pytest.raises(HTTPException) as excinfo:
        client.post("/molecules", json={"smiles": "c1ccccc1"})
    assert excinfo.value.status_code == 409
    assert excinfo.value.detail == "Molecule already exists"

#Test 4 updating a molecule by smiles
@patch.object(MoleculeDAO, 'update', new_callable=AsyncMock)
def test_update_molecule(mock_update):
    mock_update.return_value = {"id": 1, "smiles": "c1ccccc1"}
    response = client.put("/molecules/1", json={"smiles": "c1ccccc1"})
    
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "c1ccccc1"}
    mock_update.assert_called_once_with(1, MoleculeUpdate(smiles="c1ccccc1"))

#Test 5 updating a molecule by identifier that does not exist
@patch.object(MoleculeDAO, 'update', new_callable=AsyncMock)
def test_update_molecule_incorrect_id(mock_update):
    mock_update.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.put("/molecules/2", json={"smiles": "c1ccccc1"})
    assert excinfo.value.status_code == 404
    assert excinfo.value.detail == "Molecule not found"

#Test 6 updating a molecule with an invalid smiles structure
@patch.object(MoleculeDAO, 'add_molecule', new_callable=AsyncMock)
@patch.object(MoleculeDAO, 'update', new_callable=AsyncMock)
def test_update_molecule_incorrect_smiles(mock_update, mock_add_molecule):
    mock_add_molecule.return_value = 1
    client.post("/molecules", json={"smiles": "c1ccccc1"})
    
    mock_update.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.put("/molecules/1", json={"smiles": "test"})
    
    assert excinfo.value.status_code == 400
    assert excinfo.value.detail == "Invalid SMILES molecule"


#Test 7 updating a molecule with a smiles that already exists
@patch.object(MoleculeDAO, 'add_molecule', new_callable=AsyncMock)
@patch.object(MoleculeDAO, 'update', new_callable=AsyncMock)
@patch.object(MoleculeDAO, 'find_full_data', new_callable=AsyncMock)
def test_update_molecule_duplicate_smiles(mock_find_full_data, mock_update, mock_add_molecule):
    mock_add_molecule.return_value = 1
    client.post("/molecules", json={"smiles": "c1ccccc1"})
    mock_find_full_data.return_value = {"id": 1, "smiles": "c1ccccc1"}
    
    mock_update.side_effect = HTTPException(status_code=409, detail="Molecule with this SMILES value already exists")
    with pytest.raises(HTTPException) as excinfo:
        client.put("/molecules/1", json={"smiles": "c1ccccc1"})
    
    assert excinfo.value.status_code == 409
    assert excinfo.value.detail == "Molecule with this SMILES value already exists"

# Test 8 deleting a molecule
@patch.object(MoleculeDAO, 'delete_molecule_by_id', new_callable=AsyncMock)
def test_delete_molecule(mock_delete_molecule_by_id):
    mock_delete_molecule_by_id.return_value = 1
    response = client.delete("/molecules/1")
    
    assert response.status_code == 200
    assert response.json() == {"message": "The molecule with id 1 is deleted!"}
    mock_delete_molecule_by_id.assert_called_once_with(molecule_id=1)

#Test 9 deleting an identifier that doesnt exist
@patch.object(MoleculeDAO, 'delete_molecule_by_id', new_callable=AsyncMock)
def test_delete_molecule_incorrect_id(delete_molecule_by_id):
    delete_molecule_by_id.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.delete("/molecules/1")
    assert excinfo.value.status_code == 404
    assert excinfo.value.detail == "Molecule with id 1 does not exist!"

#Tests with Redis

# Test 10 retreiving molecule by identifier, not using cache

@patch.object(MoleculeDAO, 'find_full_data', new_callable=AsyncMock)
@patch("src.cache.redis_client")
@pytest.mark.asyncio
async def test_get_molecule_by_id(mock_redis_client, mock_find_full_data):
    mock_redis = mock_redis_client.return_value
    mock_redis.get = AsyncMock(return_value=None) 
    mock_redis.setex = AsyncMock()
    mock_redis.close = AsyncMock()

    mock_molecule = MagicMock()
    mock_molecule.model_dump.return_value = {"id": 1, "smiles": "c1ccccc1"}

    mock_find_full_data.return_value = mock_molecule

    response = client.get("/molecules/1")

    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "c1ccccc1"}

# Test 11 retreiving molecule by identifier, using cache

@patch.object(MoleculeDAO, 'find_full_data', new_callable=AsyncMock)
@patch("src.cache.redis_client")
@pytest.mark.asyncio
async def test_cache_usage_on_repeated_requests(mock_redis_client, mock_find_full_data):
    mock_redis = mock_redis_client.return_value
    mock_redis.get = AsyncMock(return_value=b'{"id": 1, "smiles": "c1ccccc1"}')
    mock_redis.setex = AsyncMock()
    mock_redis.close = AsyncMock()

    mock_molecule = MagicMock()
    mock_molecule.model_dump.return_value = {"id": 1, "smiles": "c1ccccc1"}

    mock_find_full_data.return_value = mock_molecule

    response = client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "c1ccccc1"}

    #Check if find_full_data was not called (cache hit) 
    mock_find_full_data.assert_not_called()

    response = client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "c1ccccc1"}

    mock_find_full_data.assert_not_called()


# Test 12 retreiving molecule by identifier that does not exist

@patch.object(MoleculeDAO, 'find_full_data', new_callable=AsyncMock)
@pytest.mark.asyncio
async def test_get_molecule_by_non_existent_id(mock_find_full_data):
    mock_find_full_data.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.get("/molecules/2")
    assert excinfo.value.status_code == 404
    assert excinfo.value.detail == "Molecule with id 2 does not exist!"

# Test 13 retreiving all molecules in DB, without caching

@patch.object(MoleculeDAO, 'find_all_molecules', new_callable=AsyncMock)
@patch("src.cache.redis_client")
@pytest.mark.asyncio
async def test_no_cache_usage_on_repeated_requests(mock_redis_client, mock_find_all_molecules):
    
    mock_redis = mock_redis_client.return_value
    mock_redis.get = AsyncMock(return_value=None)  
    mock_redis.setex = AsyncMock()
    mock_redis.close = AsyncMock()
    
    
    mock_find_all_molecules.return_value = [{"id": 1, "smiles": "c1ccccc1"}]
    
    
    response = client.get("/molecules?limit=1")
    assert response.status_code == 200
    assert response.json() == [{"id": 1, "smiles": "c1ccccc1"}]
    

# Test 14 retreiving all molecules in DB, with caching

@patch.object(MoleculeDAO, 'find_all_molecules', new_callable=AsyncMock)
@patch("src.cache.redis_client")
@pytest.mark.asyncio
async def test_cache_usage_on_retreiving_all_molecules(mock_redis_client, mock_find_all_molecules):
    mock_redis = mock_redis_client.return_value
    mock_redis.get = AsyncMock(return_value=b'{"id": 1, "smiles": "c1ccccc1"}')
    mock_redis.setex = AsyncMock()
    mock_redis.close = AsyncMock()
    
    mock_find_all_molecules.return_value = [{"id": 1, "smiles": "c1ccccc1"}]
    
    # Perform the first GET request (cache miss)
    response = client.get("/molecules?limit=1")
    assert response.status_code == 200
    assert response.json() == [{"id": 1, "smiles": "c1ccccc1"}]

    mock_find_all_molecules.assert_not_called()

    
    # Perform the second GET request (cache miss again)
    response = client.get("/molecules?limit=1")
    assert response.status_code == 200
    assert response.json() == [{"id": 1, "smiles": "c1ccccc1"}]


    mock_find_all_molecules.assert_not_called()


# Test 15 retreiving all molecules in DB, without caching

@patch.object(MoleculeDAO, 'substructure_search', new_callable=AsyncMock)
@patch("src.cache.redis_client")
@pytest.mark.asyncio
async def test_no_cache_usage_on_substructure_search(mock_redis_client, mock_substructure_search):
    mock_redis = mock_redis_client.return_value
    mock_redis.get = AsyncMock(return_value=None)  
    mock_redis.setex = AsyncMock()
    mock_redis.close = AsyncMock()
    
    mock_substructure_search.return_value = ["c1ccccc1"]
    
    response = client.get("/substructures", params={"smiles": "c1ccccc1"})
    assert response.status_code == 200
    assert response.json() == ["c1ccccc1"]
    
# Test 16 retreiving all molecules in DB, with caching

@patch.object(MoleculeDAO, 'substructure_search', new_callable=AsyncMock)
@patch("src.cache.redis_client")
@pytest.mark.asyncio
async def test_cache_usage_on_substructure_search(mock_redis_client, mock_substructure_search):
    
    mock_redis = mock_redis_client.return_value
    mock_redis.get = AsyncMock(return_value=json.dumps(["c1ccccc1"]))  
    mock_redis.setex = AsyncMock()
    mock_redis.close = AsyncMock()
    
    mock_substructure_search.return_value = ["c1ccccc1"]
    
    response = client.get("/substructures", params={"smiles": "c1ccccc1"})
    assert response.status_code == 200
    assert response.json() == ["c1ccccc1"]
    
    mock_substructure_search.assert_not_called()

    response = client.get("/substructures", params={"smiles": "c1ccccc1"})
    assert response.status_code == 200
    assert response.json() == ["c1ccccc1"]
    
    mock_substructure_search.assert_not_called()

#Test 17 substructure search with invalid smiles
@patch.object(MoleculeDAO, 'substructure_search', new_callable=AsyncMock)
@pytest.mark.asyncio
async def test_substructure_search_invalid_smiles(mock_substructure_search):
    
    mock_substructure_search.return_value = None

    with pytest.raises(HTTPException) as excinfo:
        client.put("/molecules/1", json={"smiles": "test"})
    
    assert excinfo.value.status_code == 400
    assert excinfo.value.detail == "Invalid SMILES molecule"