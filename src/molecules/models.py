from sqlalchemy.orm import Mapped
from src.database import Base, str_uniq, int_pk


class Molecule(Base):
    id: Mapped[int_pk]
    smiles: Mapped[str_uniq]

    def __str__(self):
        return (
            f"{self.__class__.__name__}(id={self.id}, "
            f"smiles={self.smiles!r})"
        )

    def __repr__(self):
        return str(self)
    
    def to_dict(self):
        return {
            'id': self.id,
            'smiles': self.smiles
        }
