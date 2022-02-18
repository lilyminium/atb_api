from typing import List, Optional, Dict
import datetime

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
    select_output_atom_names: Optional[Dict[int, str]] = None
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

    def to_mdanalysis(self, united: bool = False, optimized: bool = True):
        require_package("MDAnalysis")
        import MDAnalysis as mda
        from .mdanalysis import (
            PartialCharge,
            United,
            OutputAtomisticID,
            OutputUnitedID,
        )

        atoms = self.atoms
        if united:
            atoms = [atom for atom in atoms if atom.get_united_ljsym()]

        u = mda.Universe.empty(len(atoms), trajectory=True)
        for attr in (
            "names",
            "resnames",
            "elements",
            "types",
            "charges",
            "partial_charges",
            "united",
            "aromaticities",
            "ids",
            "output_atomistic_ids",
            "output_united_ids",
            "masses",
        ):
            u.add_TopologyAttr(attr)

        id_to_index = {}
        for i, (atbatom, mdaatom) in enumerate(zip(atoms, u.atoms)):
            mdaatom.name = atbatom.name
            mdaatom.element = atbatom.element.symbol
            mdaatom.type = atbatom.atomistic_lj_atom_type
            mdaatom.charge = atbatom.formal_charge
            mdaatom.partial_charge = atbatom.atomistic_partial_charge
            mdaatom.united = atbatom.is_united
            mdaatom.aromaticity = atbatom.is_aromatic
            mdaatom.id = i + 1
            mdaatom.output_atomistic_id = atbatom.atomistic_output_id
            mdaatom.output_united_id = atbatom.united_output_id
            mdaatom.mass = atbatom.atomistic_mass
            if optimized:
                mdaatom.position = atbatom.optimized_coordinate * 10
            else:
                mdaatom.position = atbatom.original_coordinate * 10
            mdaatom.residue.resname = atbatom.residue_name
            id_to_index[atbatom.input_id] = i

        bond_values = []
        bond_types = []
        bond_orders = []

        for bond in self.bonds:
            if not all(x in id_to_index for x in bond.atomistic_atom_ids):
                continue
            i_, j_ = bond.atomistic_atom_ids
            i = id_to_index[i_]
            j = id_to_index[j_]
            bond_values.append((i, j))
            bond_types.append((u.atoms[i].type, u.atoms[j].type))
            bond_orders.append(bond.order)

        u.add_bonds(bond_values, types=bond_types, order=bond_orders)

        for parameter_name in ("angles", "dihedrals"):  # , "impropers"):
            atb_parameters = getattr(self, parameter_name)
            values = []
            types = []
            for parameter in atb_parameters:
                if parameter_name == "dihedrals" and not parameter.essential:
                    continue
                if not all(x in id_to_index for x in parameter.atomistic_atom_ids):
                    continue
                value_ = tuple(id_to_index[x] for x in parameter.atomistic_atom_ids)
                type_ = tuple(u.atoms[x].type for x in value_)
                values.append(value_)
                types.append(type_)
            u._add_topology_objects(parameter_name, values, types=types)

        return u

    def to_itp(
        self, filename: str,
        use_input_order: bool = False,
        united: bool = False,
    ):
        itp_string = self.to_itp_string(use_input_order=use_input_order, united=united)
        with open(str(filename), "w") as f:
            f.write(itp_string)

    def to_itp_string(self, use_input_order: bool = False, united: bool = False) -> str:
        from .templates.itp import ITP_TEMPLATE

        atom_id_mapping = self.get_atom_id_mapping(
            use_input_order=use_input_order, united=united
        )
        atoms = sorted(
            [a for a in self.atoms if a.input_id in atom_id_mapping],
            key=lambda a: atom_id_mapping[a.input_id],
        )
        atom_str = "\n".join(
            [
                atom.to_itp_string(
                    output_id=atom_id_mapping[atom.input_id],
                    united=united,
                    residue_name=self.residue_name,
                )
                for atom in atoms
            ]
        )
        if united:
            charge = sum([atom.united_partial_charge for atom in atoms])
        else:
            charge = sum([atom.atomistic_partial_charge for atom in atoms])

        bonds = self.get_sorted_parameters(self.bonds, atom_id_mapping, (0, 1))
        bond_str = "\n".join([bond.to_itp_string(atom_id_mapping) for bond in bonds])
        angles = self.get_sorted_parameters(self.angles, atom_id_mapping, (1, 0, 2))
        angle_str = "\n".join(
            [angle.to_itp_string(atom_id_mapping) for angle in angles]
        )
        dihedrals = self.get_sorted_parameters(
            self.dihedrals, atom_id_mapping, (0, 1, 2, 3)
        )
        essential_dihedrals = [x for x in dihedrals if x.essential]
        dih_str = "\n".join(
            [dih.to_itp_string(atom_id_mapping) for dih in essential_dihedrals]
        )
        impropers = self.get_sorted_parameters(
            self.impropers, atom_id_mapping, (0, 1, 2, 3)
        )
        if not united:
            impropers = [
                imp
                for imp in impropers
                if not self.atoms[imp.atomistic_atom_ids[0]].is_united
            ]
        imp_str = "\n".join([imp.to_itp_string(atom_id_mapping) for imp in impropers])

        aa_output_id_to_atom = {atom.atomistic_output_id: atom for atom in atoms}
        aa_output_id_to_new_id = {
            k: v.input_id for k, v in aa_output_id_to_atom.items()
        }
        exclusions = []
        exclusions_to_include = []
        for atom1 in atoms:
            id1 = atom_id_mapping[atom1.input_id]
            excluded_atomistic_output_ids = atom1.get_exclusions(united=False)
            for atom2id in excluded_atomistic_output_ids:
                if atom2id in aa_output_id_to_new_id:
                    excl = tuple(sorted([id1, aa_output_id_to_new_id[atom2id]]))
                    exclusions.append(excl)
                    atom2 = aa_output_id_to_atom[atom2id]
                    if atom1.is_aromatic and atom2.is_aromatic:
                        for ring in self.rings:
                            if len(ring.atomistic_atom_ids) == 6:
                                if (
                                    atom1.input_id in ring.atomistic_atom_ids
                                    and atom2.input_id in ring.atomistic_atom_ids
                                ):
                                    exclusions_to_include.append(excl)

        exclusion_str = "\n".join(
            [f"{x[0]:>6d} {x[1]:>6d}" for x in sorted(set(exclusions_to_include))]
        )

        pairs = []
        for dih in dihedrals:
            first, _, __, last = dih.atomistic_atom_ids
            pair = sorted([atom_id_mapping[first], atom_id_mapping[last]])
            if pair not in exclusions:
                pairs.append(pair)
        pair_str = "\n".join([f"{x[0]:>5d}{x[1]:>5d}    1" for x in sorted(pairs)])

        now = datetime.datetime.now()
        resolution = "all atom" if not united else "united atom"

        return ITP_TEMPLATE.format(
            time=now.strftime("%H:%M"),
            date=now.strftime("%Y-%m-%d"),
            revision=self.REV_DATE,
            residue_name=self.residue_name,
            resolution_upper=resolution.upper(),
            resolution=resolution,
            molecule_molid=self.molid,
            molecule_hash=self.topology_hash,
            atoms=atom_str,
            total_charge=charge,
            bonds=bond_str,
            angles=angle_str,
            dihedrals=dih_str,
            impropers=imp_str,
            exclusions=exclusion_str,
            pairs=pair_str,
        )

    def get_sorted_parameters(self, initial_container, atom_id_mapping, sort_indices):
        return sorted(
            [
                b
                for b in initial_container
                if all(a in atom_id_mapping for a in b.atomistic_atom_ids)
            ],
            key=lambda b: tuple(
                atom_id_mapping[b.atomistic_atom_ids[i]] for i in sort_indices
            ),
        )

    def get_atom_id_mapping(
        self, use_input_order: bool = False, united: bool = False
    ) -> Dict[int, int]:
        id_to_atoms = {atom.input_id: atom for atom in self.atoms}

        if use_input_order:
            if united:
                id_to_atoms = {k: v for k, v in id_to_atoms.items() if v.get_united_ljsym()}
            id_to_output_id = {k: i for i, k in enumerate(sorted(id_to_atoms), 1)}
        else:
            if united:
                id_to_output_id = {
                    k: v.united_output_id for k, v in id_to_atoms.items()
                    if v.get_united_ljsym()
                }
            else:
                id_to_output_id = {
                    k: v.atomistic_output_id for k, v in id_to_atoms.items()
                }
        return id_to_output_id
