from fastapi import HTTPException, status
from sqlalchemy import delete, update
from typing import Optional 
from sqlalchemy.future import select
from sqlalchemy.exc import NoResultFound, IntegrityError
from src.molecules.models import Molecule
from src.database import async_session_maker
from src.dao.base import BaseDAO
from src.molecules.schema import MoleculeResponse, MoleculeUpdate
from rdkit import Chem



class MoleculeDAO(BaseDAO):
    model = Molecule

    @classmethod
    async def find_all_molecules(cls, limit: Optional[int] = None):
        async with async_session_maker() as session:
            query = select(cls.model)
            if limit:
                query = query.limit(limit)
            result = await session.execute(query)
            molecules = result.scalars().all()
            return [molecule.to_dict() for molecule in molecules]

    @classmethod
    async def add_molecule(cls, **molecule_data: dict):
        async with async_session_maker() as session:
            try:
                async with session.begin():
                    new_molecule = cls.model(**molecule_data)
                    session.add(new_molecule)
                    await session.flush()
                    new_molecule_id = new_molecule.id
                    await session.commit()
                    return new_molecule_id
            except IntegrityError as e:
                if 'molecules_smiles_key' in str(e):
                    raise HTTPException(
                        status_code=status.HTTP_409_CONFLICT,
                        detail="Molecule already exists"
                    )
                else:
                    raise e
            
    @classmethod
    async def find_full_data(cls, molecule_id: int):
        async with async_session_maker() as session:
            query = select(cls.model).filter_by(id=molecule_id)
            result = await session.execute(query)
            molecules_info = result.scalar_one_or_none()

            if not molecules_info:
                return None

                       
            return MoleculeResponse(id=molecules_info.id, smiles=molecules_info.smiles)



    @classmethod
    async def delete_molecule_by_id(cls, molecule_id: int):
        async with async_session_maker() as session:
            async with session.begin():
                query = select(cls.model).filter_by(id=molecule_id)
                result = await session.execute(query)
                molecule_to_delete = result.scalar_one_or_none()

                if not molecule_to_delete:
                    return None

                await session.execute(delete(cls.model).filter_by(id=molecule_id))

                await session.commit()
                return molecule_id
 
    @classmethod
    async def update(cls, molecule_id: int, update_data: MoleculeUpdate) -> MoleculeResponse:
        async with async_session_maker() as session:
            async with session.begin():
                stmt = select(cls.model).where(cls.model.id == molecule_id)
                result = await session.execute(stmt)
                
                try:
                    molecule = result.scalar_one()
                except NoResultFound:
                    return None
                
                try:
                    stmt = (
                        update(cls.model)
                        .where(cls.model.id == molecule_id)
                        .values(smiles=update_data.smiles)
                        .returning(cls.model)
                    )
                    result = await session.execute(stmt)
                    updated_molecule = result.scalar_one()
                except IntegrityError as e:
                    if 'molecules_smiles_key' in str(e):
                        raise HTTPException(
                            status_code=status.HTTP_409_CONFLICT,
                            detail="Molecule with this SMILES value already exists"
                        )
                    else:
                        raise HTTPException(
                            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail="Database error"
                        )
                
                return MoleculeResponse(id=updated_molecule.id, smiles=updated_molecule.smiles)
    
    @classmethod
    async def substructure_search(cls, smiles: str):
        substructure = Chem.MolFromSmiles(smiles)
        if not substructure:
            raise ValueError("Invalid SMILES molecule")

        async with async_session_maker() as session:
            query = select(cls.model)
            results = await session.execute(query)
            molecules = results.scalars().all()

            matches = []
            for molecule in molecules:
                mol = Chem.MolFromSmiles(molecule.smiles)
                if mol and mol.HasSubstructMatch(substructure):
                    matches.append(molecule.smiles)

            return matches