class RBMolecule:
    def __init__(
        self,
        molecule_id: int | None = None,
        smiles: str | None = None,
    ):
        self.id = molecule_id
        self.smiles = smiles

    def to_dict(self) -> dict:
        data = {
            "id": self.id,
            "smiles": self.smiles,
        }
        filtered_data = {key: value for key, value in data.items() if value is not None}
        return filtered_data

