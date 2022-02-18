from typing import List, Dict, Any, Optional

import numpy as np

from elementable import Elementable
from pydantic import Field, validator

from .base import Model


class ElementBase(Model):

    pass


elements = Elementable(element_cls=ElementBase)
Element = elements.element_class


class Atom(Model):
    """ATB Atom"""

    is_aromatic: bool = Field(default=False, alias="aromatic")
    is_united: bool = Field(default=False, alias="united")

    charge_group: int = Field(default=0, alias="cgroup")
    equivalence_group: int = Field(default=0, alias="equivalenceGroup")
    formal_charge: int = Field(default=0)
    residue_name: str = Field(default="UNK", alias="group")
    standard_atom_type_code: int = Field(default=0, alias="std_iacm")
    input_id: int = Field(
        default=1,
        alias="id",
        description="Index + 1 for ordered index of atoms in the input PDB",
    )
    name: str = Field(default="", alias="symbol")
    element: Element = Field(default=elements.X, alias="type")
    type_energy: str = ""

    atomistic_partial_charge: float = Field(default=0, alias="charge")
    atomistic_bonded_atoms: List[int] = Field(default_factory=list, alias="conn")
    atomistic_exclusion_atoms: List[int] = Field(default_factory=list, alias="excl")
    atomistic_atb_atom_type_code: int = Field(default=1, alias="iacm")
    atomistic_charge_group_code: int = Field(default=1, alias="icgm")
    atomistic_output_id: int = Field(
        default=1,
        alias="index",
        description="Index + 1 for ordered index of atoms in the output files",
    )
    atomistic_lj_atom_type: str = Field(
        default="", alias="ljsym", description="Nonbonded Atom type for LJ parameters"
    )
    atomistic_mass: float = Field(default=0, alias="mass")
    atomistic_mass_code: int = Field(default=0, alias="mass_code")

    united_partial_charge: float = Field(default=None, alias="ucharge")
    united_bonded_atoms: List[int] = Field(default_factory=list, alias="uconn")
    united_exclusion_atoms: List[int] = Field(default_factory=list, alias="uexcl")
    united_atb_atom_type_code: int = Field(default=1, alias="uiacm")
    united_charge_group_code: int = Field(default=1, alias="uicgm")
    united_output_id: int = Field(
        default=1,
        alias="uindex",
        description="Index + 1 for ordered index of atoms in the output files",
    )
    united_lj_atom_type: str = Field(
        default=None, alias="uljsym", description="Nonbonded Atom type for LJ parameters"
    )
    united_mass: float = Field(default=None, alias="umass")
    united_mass_code: int = Field(default=0, alias="umass_code")

    original_coordinate: np.ndarray = Field(
        default=[0, 0, 0], alias="coord", description="coordinates in nm"
    )
    optimized_coordinate: np.ndarray = Field(
        default=[0, 0, 0], alias="ocoord", description="coordinates in nm"
    )

    pdb: str = ""

    @validator("element", pre=True)
    def _validate_element_from_symbol(cls, v):
        if isinstance(v, str):
            v = elements(symbol=v.capitalize())
        return v

    @validator("original_coordinate", "optimized_coordinate", pre=True)
    def _validate_coordinates(cls, v):
        v = np.asarray(v)
        return v

    @validator("united_partial_charge", "united_mass", "united_lj_atom_type", pre=True)
    def _set_to_atomistic(cls, v, values, field):
        if v is None and ("united" in values or "is_united" in values):
            united = values.get("is_united", values.get("united"))
            if not united:
                field_name = field.name.replace("united", "atomistic")
                at_field = cls.__fields__[field_name]
                if at_field.alias in values:
                    v = values[at_field.alias]
                if field_name in values:
                    v = values[field_name]
            elif field.type_ == str:
                v = ""
            else:
                v = 0
        return v

    @classmethod
    def from_atb_dict(cls, atb_dict: Dict[str, Any] = {}):
        return cls(**atb_dict)

    def to_itp_string(
        self,
        output_id: Optional[int] = None,
        united: Optional[bool] = None,
        residue_name: Optional[str] = None,
    ):
        if residue_name is None:
            residue_name = self.residue_name

        if united is None:
            united = self.is_united

        if not united:
            if output_id is None:
                output_id = self.atomistic_output_id
            ljtype = self.atomistic_lj_atom_type
            charge = self.atomistic_partial_charge
            mass = self.atomistic_mass

        else:
            if output_id is None:
                output_id = self.united_output_id
            ljtype = self.get_united_ljsym()
            charge = self.united_partial_charge
            mass = self.united_mass

        return (
            f"{output_id:>5d} {ljtype:>5}    1 "
            f"{residue_name:>7} {self.name:>6} {output_id:>4d} "
            f"{charge:>8.3f} {mass:>8.4f}"
        )

    def get_exclusions(self, united: Optional[bool] = None):
        if united is None:
            united = self.is_united

        if united:
            return self.united_exclusion_atoms
        return self.atomistic_exclusion_atoms

    def get_united_ljsym(self):
        if not self.united_lj_atom_type and not self.is_united:
            return self.atomistic_lj_atom_type
        return self.united_lj_atom_type
