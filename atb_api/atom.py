from typing import List, Dict, Any

import numpy as np

from elementable import Element, Elements
from pydantic import Field, validator

from .base import Model


class Atom(Model):
    """ATB Atom"""
    is_aromatic: bool = Field(default=False, alias="aromatic")
    is_united: bool = Field(default=False, alias="united")
    charge_group: int = Field(default=0, alias="cgroup")
    equivalence_group: int = Field(default=0, alias="equivalenceGroup")
    formal_charge: int = Field(default=0)
    residue_name: str = Field(default="UNK", alias="group")
    standard_atom_type_code: int = Field(default=0, alias="std_iacm")
    input_id: int = Field(default=1, alias="id", description="Index + 1 for ordered index of atoms in the input PDB")
    name: str = Field(default="", alias="symbol")
    element: Element = Field(default=Elements.X, alias="type")
    type_energy: str = ""

    atomistic_partial_charge: float = Field(default=0, alias="charge")
    atomistic_bonded_atoms: List[int] = Field(default_factory=list, alias="conn")
    atomistic_exclusion_atoms: List[int] = Field(default_factory=list, alias="excl")
    atomistic_atb_atom_type_code: int = Field(default=1, alias="iacm")
    atomistic_charge_group_code: int = Field(default=1, alias="icgm")
    atomistic_output_id: int = Field(default=1, alias="index",
                                     description="Index + 1 for ordered index of atoms in the output files")
    atomistic_lj_atom_type: str = Field(default="", alias="ljsym", description="Nonbonded Atom type for LJ parameters")
    atomistic_mass: float = Field(default=0, alias="mass")
    atomistic_mass_code: int = Field(default=0, alias="mass_code")

    united_partial_charge: float = Field(default=0, alias="ucharge")
    united_bonded_atoms: List[int] = Field(default_factory=list, alias="uconn")
    united_exclusion_atoms: List[int] = Field(default_factory=list, alias="uexcl")
    united_atb_atom_type_code: int = Field(default=1, alias="uiacm")
    united_charge_group_code: int = Field(default=1, alias="uicgm")
    united_output_id: int = Field(default=1, alias="uindex",
                                  description="Index + 1 for ordered index of atoms in the output files")
    united_lj_atom_type: str = Field(default="", alias="uljsym", description="Nonbonded Atom type for LJ parameters")
    united_mass: float = Field(default=0, alias="umass")
    united_mass_code: int = Field(default=0, alias="umass_code")

    original_coordinate: np.ndarray = Field(default=[0, 0, 0], alias="coord", description="coordinates in nm")
    optimized_coordinate: np.ndarray = Field(default=[0, 0, 0], alias="ocoord", description="coordinates in nm")

    pdb: str = ""

    @validator("element", pre=True)
    def _validate_element_from_symbol(cls, v):
        if isinstance(v, str):
            v = Element(symbol=v.capitalize())
        return v

    @validator("original_coordinate", "optimized_coordinate", pre=True)
    def _validate_coordinates(cls, v):
        v = np.asarray(v)
        return v

    @classmethod
    def from_atb_dict(cls, atb_dict: Dict[str, Any] = {}):
        return cls(**atb_dict)
