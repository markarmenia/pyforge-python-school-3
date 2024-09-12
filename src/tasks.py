import asyncio
from src.celery_worker import celery_app 
from src.molecules.dao import MoleculeDAO
from rdkit import Chem

@celery_app.task(name="substructure_search_task")
def substructure_search_task(smiles: str):
    substructure = Chem.MolFromSmiles(smiles)
    if not substructure:
        raise ValueError("Invalid SMILES structure")
    
    # Run the asynchronous substructure search synchronously
    matches = asyncio.run(MoleculeDAO.substructure_search(smiles))
    
    if not matches:
        return {"message": "No results found"}
    
    return matches
