from typing import List, Optional, Dict

import numpy as np
from pydantic import Field, validator

from .base import Model
from .atom import Atom
from .parameter import Bond, Angle, Dihedral, Improper, Ring
from .utils import require_package


class Molecule(Model):
    """ATB Molecule"""
    molid: int
    name: Optional[str] = None
    iupac: Optional[str] = None
    residue_name: str = Field(alias="rnme")
    topology_hash: str = ""
    total_charge: int = 0
    total_qm_charge: float = 0
    update_ifp: bool = False
    use_CH1_united_atom_for_double_bonds: bool = False
    use_charge_assign_charge_model: bool = True
    run_charge_group_partitioning: bool = False
    enforce_unique_atom_names: bool = True
    dipole: np.ndarray
    scale_manual_charges: bool = True
    symmetrize_charges: bool = Field(default=True, alias="symmetrise_charges")
    sp_timeout: Optional[int] = None
    select_output_atom_names: Dict[int, str] = {}
    ff_version: str
    charge_assign_error: str
    atom_order_template: Optional[str] = None
    amino_acid_building_block: bool = False
    new_parameter_added: bool = False

    additional_info_lines: List[str] = []
    VERSION: str
    REV_DATE: str
    A_born_energy: float
    E_born_energy: float

    qm_level: int = 1

    qm0_method: str
    qm_1_ESP: np.ndarray
    manual_charges: Optional[List[float]] = None
    has_manual_charges: bool = False

    volume: float

    atoms: List[Atom] = []
    bonds: List[Bond] = []
    angles: List[Angle] = []
    dihedrals: List[Dihedral] = []
    impropers: List[Improper] = []
    rings: List[Ring] = []

    @validator("qm_1_ESP", "dipole", pre=True)
    def _to_numpy_array(cls, v):
        v = np.asarray(v)
        return v

    @validator("atoms", "rings", pre=True)
    def _to_list(cls, v):
        if isinstance(v, dict):
            v = [v[x] for x in sorted(v)]
        return v

    @classmethod
    def from_atb_dict(cls, v):
        v = dict(v)
        dct = v.pop("var", {})
        v.pop("_dihedrals", None)
        dct.update(v)
        return cls(**dct)

    def to_rdkit(self):
        require_package("rdkit")
        from rdkit import Chem

        rdmol = Chem.RWMol()
        atbatoms = sorted(self.atoms, key=lambda x: x.input_id)
        rdconf = Chem.Conformer(len(atbatoms))
        optimized_coordinates = []
        id_to_index = {}
        for i, atbatom in enumerate(atbatoms):
            rdatom = Chem.Atom(atbatom.element.symbol)
            rdatom.SetIsAromatic(atbatom.is_aromatic)
            rdatom.SetFormalCharge(atbatom.formal_charge)
            rdatom.SetAtomMapNum(atbatom.atomistic_output_id)

            # coord: original, nm
            original = atbatom.original_coordinate * 10
            rdatom.SetDoubleProp("original_x", original[0])
            rdatom.SetDoubleProp("original_y", original[1])
            rdatom.SetDoubleProp("original_z", original[2])
            # ocoord: optimized, nm
            optimized = atbatom.optimized_coordinate * 10
            rdatom.SetDoubleProp("optimized_x", optimized[0])
            rdatom.SetDoubleProp("optimized_y", optimized[1])
            rdatom.SetDoubleProp("optimized_z", optimized[2])

            rdatom.SetBoolProp("is_united", atbatom.is_united)

            optimized_coordinates.append(optimized)

            rdconf.SetAtomPosition(i, original)
            rdmol.AddAtom(rdatom)
            id_to_index[atbatom.input_id] = i

        BONDTYPES = {
            1.0: Chem.BondType.SINGLE,
            1.5: Chem.BondType.AROMATIC,
            2.0: Chem.BondType.DOUBLE,
            3.0: Chem.BondType.TRIPLE,
        }

        for bond in self.bonds:
            i = id_to_index[bond.atomistic_atom_ids[0]]
            j = id_to_index[bond.atomistic_atom_ids[1]]
            index = rdmol.AddBond(i, j, BONDTYPES[bond.order]) - 1
            rdbond = rdmol.GetBondWithIdx(index)
            rdbond.SetIsAromatic(bond.is_aromatic)
            rdbond.SetBoolProp("is_united", bond.is_united)

        Chem.SanitizeMol(rdmol)
        return Chem.Mol(rdmol)
