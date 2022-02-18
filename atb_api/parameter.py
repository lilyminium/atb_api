from typing import List, Optional, Dict, Tuple, ClassVar

from pydantic import Field, validator

from .base import Model


class ParameterCode(Model):
    value: float
    code: int
    closest: bool = False
    nostd: bool = False
    match_r: List[str] = []
    fc: float
    types: List[List[str]] = []
    uncertain: bool = False

    @validator("types", pre=True)
    def _validate_types(cls, v):
        if len(v) and not isinstance(v[0], list):
            v = [v]
        return v


class Parameter(Model):

    atomistic_atom_ids: List[int] = Field(alias="atoms")
    is_united: bool = Field(default=False, alias="united")
    value: float
    fc: Optional[float] = None
    hfc: Optional[float] = None
    fc_QM_hessian: Optional[float] = None
    hfc_QM_hessian: Optional[float] = None

    def atom_ids_to_string(self, atom_id_mapping: Dict[int, int]) -> str:
        atom_ids = [atom_id_mapping[x] for x in self.atomistic_atom_ids]
        return " " + " ".join([f"{x:>4d}" for x in atom_ids])

    def to_itp_string(self, atom_id_mapping: Dict[int, int]):
        if not all(aid in atom_id_mapping for aid in self.atomistic_atom_ids):
            return
        atom_ids = self.atom_ids_to_string(atom_id_mapping=atom_id_mapping)
        code = self.code[0].to_itp_string()
        return f"{atom_ids} {code}"


class BondCode(ParameterCode):
    hfc: float
    funct: ClassVar[int] = 2

    def to_itp_string(self):
        return f"{self.funct:>4d} {self.value:8.4f}   {self.fc:.4e}"


class Bond(Parameter):
    is_aromatic: bool = Field(default=False, alias="aromatic")
    order: float = Field(default=1)
    order_qm: float = Field(default=1)
    code: List[BondCode] = []


class AngleCode(ParameterCode):
    hfc: float
    funct: ClassVar[int] = 2

    def to_itp_string(self):
        return f"{self.funct:>4d} {self.value:>9.2f} {self.fc:>8.2f}"


class Angle(Parameter):
    code: List[AngleCode] = []


class DihedralCode(ParameterCode):
    mul: int = 1
    fc: float
    funct: ClassVar[int] = 1

    def to_itp_string(self):
        return f"{self.funct:>4d} {self.value:9.2f} {self.fc:>8.2f} {self.mul:>4d}"


class Dihedral(Parameter):
    essential: bool = False
    mul: int
    code: List[DihedralCode]


class Improper(Parameter):
    code: int

    def to_itp_string(self, atom_id_mapping: Dict[int, int]):
        if not all(aid in atom_id_mapping for aid in self.atomistic_atom_ids):
            return
        atom_ids = self.atom_ids_to_string(atom_id_mapping=atom_id_mapping)
        if self.code == 1:
            code = "   2      0.00   167.36"
        elif self.code == 2:
            code = "   2     35.26   334.72"
        else:
            raise NotImplementedError("codes that aren't 1 or 2 are not supported yet")
        return f"{atom_ids} {code}"


class Ring(Model):
    atomistic_atom_ids: List[int] = Field(alias="atoms")
    aromatic: bool = False
