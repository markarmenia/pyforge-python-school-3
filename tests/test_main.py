import pytest
from fastapi import HTTPException
from fastapi.testclient import TestClient
from unittest.mock import patch, AsyncMock
from src.molecules.router import router
from src.molecules.dao import MoleculeDAO
from src.molecules.schema import MoleculeAdd, MoleculeUpdate

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

# Test 3: Adding molecule with a SMILES that already exists
@patch.object(MoleculeDAO, 'add_molecule', new_callable=AsyncMock)
def test_add_existing_molecule(mock_add_molecule):
    mock_add_molecule.side_effect = [None, HTTPException(status_code=409, detail="Molecule already exists")]
    response = client.post("/molecules", json={"smiles": "c1ccccc1"})
    assert response.status_code == 201
    with pytest.raises(HTTPException) as excinfo:
        client.post("/molecules", json={"smiles": "c1ccccc1"})
    assert excinfo.value.status_code == 409
    assert excinfo.value.detail == "Molecule already exists"

#Test 4 retreiving molecule by identifier
@patch.object(MoleculeDAO, 'find_full_data', new_callable=AsyncMock)
def test_get_molecule_by_id(mock_find_full_data):
    mock_find_full_data.return_value = {"id": 1, "smiles": "c1ccccc1"}
    response = client.get("/molecules/1")
    
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "c1ccccc1"}
    mock_find_full_data.assert_called_once_with(molecule_id=1)

#Test 5 retreiving molecule by identifier that does not exist
@patch.object(MoleculeDAO, 'find_full_data', new_callable=AsyncMock)
def test_get_molecule_by_non_existent_id(mock_find_full_data):
    mock_find_full_data.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.get("/molecules/2")
    assert excinfo.value.status_code == 404
    assert excinfo.value.detail == "Molecule with id 2 does not exist!"

#Test 6 updating a molecule by smiles
@patch.object(MoleculeDAO, 'update', new_callable=AsyncMock)
def test_update_molecule(mock_update):
    mock_update.return_value = {"id": 1, "smiles": "c1ccccc1"}
    response = client.put("/molecules/1", json={"smiles": "c1ccccc1"})
    
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "c1ccccc1"}
    mock_update.assert_called_once_with(1, MoleculeUpdate(smiles="c1ccccc1"))

#Test 7 updating a molecule by identifier that does not exist
@patch.object(MoleculeDAO, 'update', new_callable=AsyncMock)
def test_update_molecule_incorrect_id(mock_update):
    mock_update.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.put("/molecules/2", json={"smiles": "c1ccccc1"})
    assert excinfo.value.status_code == 404
    assert excinfo.value.detail == "Molecule not found"

#Test 8 updating a molecule with an invalid smiles structure
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


#Test 9 updating a molecule with a smiles that already exists
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


# Test 10 deleting a molecule
@patch.object(MoleculeDAO, 'delete_molecule_by_id', new_callable=AsyncMock)
def test_delete_molecule(mock_delete_molecule_by_id):
    mock_delete_molecule_by_id.return_value = 1
    response = client.delete("/molecules/1")
    
    assert response.status_code == 200
    assert response.json() == {"message": "The molecule with id 1 is deleted!"}
    mock_delete_molecule_by_id.assert_called_once_with(molecule_id=1)

#Test 11 deleting an identifier that doesnt exist
@patch.object(MoleculeDAO, 'delete_molecule_by_id', new_callable=AsyncMock)
def test_delete_molecule_incorrect_id(delete_molecule_by_id):
    delete_molecule_by_id.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.delete("/molecules/1")
    assert excinfo.value.status_code == 404
    assert excinfo.value.detail == "Molecule with id 1 does not exist!"

#Test 12 retrieving all molecules
@patch.object(MoleculeDAO, 'find_all_molecules', new_callable=AsyncMock)
def test_get_all_molecules(mock_find_all_molecules):
    mock_find_all_molecules.return_value = [
        {'id': 1, 'smiles': 'c1ccccc1'},
        {'id': 2, 'smiles': 'CC(=O)O'}
    ]
    response = client.get("/molecules")
    assert response.status_code == 200
    assert response.json() == [
        {'id': 1, 'smiles': 'c1ccccc1'},
        {'id': 2, 'smiles': 'CC(=O)O'}
    ]

#Test 13 substructure search
@patch.object(MoleculeDAO, 'substructure_search', new_callable=AsyncMock)
def test_substructure_search(mock_substructure_search):
    mock_substructure_search.return_value = ["c1ccccc1"]
    response = client.get("/substructures", params={"smiles": "c1ccccc1"})
    
    assert response.status_code == 200
    assert response.json() == ["c1ccccc1"]
    mock_substructure_search.assert_called_once_with("c1ccccc1")

#Test 14 substructure search with no matching results
@patch.object(MoleculeDAO, 'substructure_search', new_callable=AsyncMock)
def test_substructure_search_no_matches(mock_substructure_search):
    mock_substructure_search.return_value = None
    with pytest.raises(HTTPException) as excinfo:
        client.get("/substructures", params={"smiles": "c1ccccc1"})
    
    assert excinfo.value.status_code == 404
    assert excinfo.value.detail == "No matching molecules found"

