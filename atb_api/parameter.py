from typing import List, Optional

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


class BondCode(ParameterCode):
    hfc: float


class Bond(Parameter):
    is_aromatic: bool = Field(default=False, alias="aromatic")
    order: float = Field(default=1)
    order_qm: float = Field(default=1)
    code: List[BondCode] = []


class AngleCode(ParameterCode):
    hfc: float


class Angle(Parameter):
    code: List[AngleCode] = []


class DihedralCode(ParameterCode):
    mul: int = 1
    fc: float


class Dihedral(Parameter):
    essential: bool = False
    mul: int
    code: List[DihedralCode]


class Improper(Parameter):
    code: int


class Ring(Model):
    atomistic_atom_ids: List[int] = Field(alias="atoms")
    aromatic: bool = False
