from fastapi.testclient import TestClient
from .main import app, database
import pytest

client = TestClient(app)


def setup_function():
    database.clear()

@pytest.fixture(autouse=True)
def run_around_tests():
    setup_function()

#Test 1 adding molecule (smiles) and its identifier
@pytest.mark.parametrize("molecule_id, smiles", [
    (0, "CCO"),
    (1, "c1ccccc1"),
    (2, "CC(=O)O"),
    (3, "CC(=O)Oc1ccccc1C(=O)O")
])
def test_add_molecules(molecule_id, smiles):
    response = client.post("/molecules", json={"id": molecule_id, "smiles": smiles})
    assert response.status_code == 201
    assert response.json() == {"id": molecule_id, "name": smiles}

#Test 2 adding molecule with an identifier that already exists
@pytest.mark.parametrize("molecule_id, smiles", [
    (0, "CCO"),
    (1, "c1ccccc1"),
    (2, "CC(=O)O"),
    (3, "CC(=O)Oc1ccccc1C(=O)O")
])
def test_add_molecules_by_same_id(molecule_id, smiles):
    client.post("/molecules", json={"id": molecule_id, "smiles": smiles})
    response = client.post("/molecules", json={"id": molecule_id, "smiles": smiles})
    assert response.status_code == 409
    assert response.json() == {"detail": "Molecule identifier already exists"}

#Test 3 adding molecule with an invalid smiles structure
@pytest.mark.parametrize("molecule_id, invalid_smiles", [
    (0, "invalid"),
    (1, "test"),
    (2, "1"),
    (3, "[]")
])
def test_add_molecules_by_invalid_smiles(molecule_id, invalid_smiles):
    response = client.post("/molecules", json={"id": molecule_id, "smiles": invalid_smiles})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES molecule"}

#Test 4 retreiving molecule by identifier
@pytest.mark.parametrize("molecule_id, smiles", [
    (0, "CCO"),
    (1, "c1ccccc1"),
    (2, "CC(=O)O"),
    (3, "CC(=O)Oc1ccccc1C(=O)O")
])
def test_get_molecules_by_id(molecule_id, smiles):
    client.post("/molecules", json={'id': molecule_id, "smiles": smiles})
    response = client.get(f"/molecules/{molecule_id}")
    assert response.status_code == 200
    assert response.json() == {'id': molecule_id, 'name': smiles}


#Test 5 updating a molecule by identifier
@pytest.mark.parametrize("updated_smiles", [
    "CCO",
    "c1ccccc1",
    "CC(=O)O",
    "CC(=O)Oc1ccccc1C(=O)O"
])
def test_update_molecules(updated_smiles):
    client.post("/molecules", json={'id': 0, "smiles": "CCO"})
    response = client.put("/molecules/0", json={'id': 0, "smiles": updated_smiles})
    assert response.status_code == 200
    assert response.json() == {'id': 0, 'name':updated_smiles}

#Test 6 updating a molecule by identifier that does not exist
@pytest.mark.parametrize("incorrect_id", [
    1,2,3,4
])
def test_update_molecules_by_incorrect_id(incorrect_id):
    client.post("/molecules", json={'id': 0, "smiles": "CCO"})
    response = client.put(f"/molecules/{incorrect_id}", json={'id': incorrect_id, "smiles": "c1ccccc1"})
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}

#Test 7 updating a molecule with an invalid smiles structure
@pytest.mark.parametrize("invalid_smiles", [
    "invalid",
    "test",
    "1",
    "[]"
])
def test_update_molecules_by_incorrect_smiles(invalid_smiles):
    client.post("/molecules", json={'id': 0, "smiles": "CCO"})
    response = client.put("/molecules/0", json={'id': 0, "smiles": invalid_smiles})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES molecule"}

#Test 8 deleting a molecule
@pytest.mark.parametrize("molecule_id, smiles", [
    (0, "CCO"),
    (1, "c1ccccc1"),
    (2, "CC(=O)O"),
    (3, "CC(=O)Oc1ccccc1C(=O)O")
])
def test_delete_molecules(molecule_id, smiles):
    client.post("/molecules", json={'id': molecule_id, "smiles": smiles})
    response = client.delete(f"/molecules/{molecule_id}")
    assert response.status_code == 200
    assert response.json() == {'id': molecule_id, 'name':smiles}


#Test 9 deleting a molecule with identifier that does not exist
@pytest.mark.parametrize("incorrect_id", [
    1,2,3,4
])
def test_delete_molecules_by_invalid_identifier(incorrect_id):
    client.post("/molecules", json={'id': 0, "smiles": "CCO"})
    response = client.delete(f"/molecules/{incorrect_id}")
    assert response.status_code == 404
    assert response.json() == {"detail":"Molecule not found"}

#Test 10 retrieving all molecules

@pytest.mark.parametrize("molecules, expected_result", [
    ([{'id': 0, "smiles": "CCO"}],
     [{'id': 0, 'name': "CCO"}]),
    
    ([{'id': 0, "smiles": "CCO"}, {'id': 1, "smiles": "c1ccccc1"}],
     [{'id': 0, 'name': "CCO"}, {'id': 1, "name": "c1ccccc1"}]),
    
    ([{'id': 0, "smiles": "CCO"}, {'id': 1, "smiles": "c1ccccc1"}, {'id': 2, "smiles": "CC(=O)O"}],
     [{'id': 0, 'name': "CCO"}, {'id': 1, "name": "c1ccccc1"}, {'id': 2, "name": "CC(=O)O"}]),
    
    ([{'id': 0, "smiles": "CCO"}, {'id': 1, "smiles": "c1ccccc1"}, {'id': 2, "smiles": "CC(=O)O"},
      {'id': 3, "smiles": "CC(=O)Oc1ccccc1C(=O)O"}],
     [{'id': 0, 'name': "CCO"}, {'id': 1, "name": "c1ccccc1"}, {'id': 2, "name": "CC(=O)O"},
      {'id': 3, "name": "CC(=O)Oc1ccccc1C(=O)O"}])
])
def test_get_molecules(molecules, expected_result):
    for molecule in molecules:
        response = client.post("/molecules", json=molecule)
        assert response.status_code == 201
    response = client.get("/molecules")
    assert response.status_code == 200
    assert response.json() == expected_result


#Test 11 substructure search
@pytest.mark.parametrize("molecules, smiles, result", [
    ([{'id': 0, "smiles": "CCO"}, {'id': 1, "smiles": "c1ccccc1"}, {'id': 2, "smiles": "CC(=O)O"},
      {'id': 3, "smiles": "CC(=O)Oc1ccccc1C(=O)O"}], "CCO", ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),
    ([{'id': 0, "smiles": "CCO"}, {'id': 1, "smiles": "c1ccccc1"}, {'id': 2, "smiles": "CC(=O)O"},
      {'id': 3, "smiles": "CC(=O)Oc1ccccc1C(=O)O"}], "c1ccccc1", ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),
    ([{'id': 0, "smiles": "CCO"}, {'id': 1, "smiles": "c1ccccc1"}, {'id': 2, "smiles": "CC(=O)O"},
      {'id': 3, "smiles": "CC(=O)Oc1ccccc1C(=O)O"}], "CC(=O)O", ["CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),
    ([{'id': 0, "smiles": "CCO"}, {'id': 1, "smiles": "c1ccccc1"}, {'id': 2, "smiles": "CC(=O)O"},
      {'id': 3, "smiles": "CC(=O)Oc1ccccc1C(=O)O"}], "CC(=O)Oc1ccccc1C(=O)O", ["CC(=O)Oc1ccccc1C(=O)O"])
])
def test_substructure_search(molecules, smiles, result):
    for molecule in molecules:
        response = client.post("/molecules", json=molecule)
        assert response.status_code == 201
    response = client.get("/substructures", params={"molecule": smiles})
    assert response.status_code == 200
    assert response.json() == result

# #Test 12 substructure search with invalid smiles
@pytest.mark.parametrize("invalid_smiles", [
    "invalid",
    "test",
    "1",
    "[]"
])
def test_substructure_search_invalid_smiles(invalid_smiles):
    client.post("/molecules", json={'id': 0, "smiles": "CCO"})
    response = client.get("/substructures", params={"molecule": invalid_smiles})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES molecule"}