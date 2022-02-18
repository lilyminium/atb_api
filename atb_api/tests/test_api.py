import pytest

from atb_api.api import ATBApi
from atb_api.exceptions import InvalidTokenError
from .datafiles import (
    MOL_903922_ITP_AA,
    MOL_903922_ITP_UA,
    MOL_903922_PDB_OPT_AA,
    MOL_903922_PDB_OPT_UA,
    MOL_903922_PDB_ORIGINAL_UA,
    MOL_903922_PDB_ORIGINAL_AA,
    MOL_903922_ITP_AA_INPUT_ORDER,
    MOL_903922_ITP_UA_INPUT_ORDER,
)


def test_invalid_authentication():
    with pytest.raises(InvalidTokenError):
        ATBApi(api_token="invalid token")


def test_submit_existing_molecule(api):
    molid = api.submit_molecule(MOL_903922_PDB_ORIGINAL_AA)
    assert molid == 903922


@pytest.mark.parametrize(
    "format, resolution, optimized, reference",
    [
        ("itp", "united_atom", None, MOL_903922_ITP_UA),
        ("itp", "ua", None, MOL_903922_ITP_UA),
        ("itp", "all_atom", None, MOL_903922_ITP_AA),
        ("itp", "aa", None, MOL_903922_ITP_AA),
        ("pdb", "united_atom", True, MOL_903922_PDB_OPT_UA),
        ("pdb", "all_atom", True, MOL_903922_PDB_OPT_AA),
        ("pdb", "united_atom", False, MOL_903922_PDB_ORIGINAL_UA),
        ("pdb", "all_atom", False, MOL_903922_PDB_ORIGINAL_AA),
    ],
)
def test_download_molecule_file(api, format, resolution, optimized, reference):
    output = api.download_molecule_file(
        molid=903922, format=format, resolution=resolution, optimized=optimized
    )
    with open(reference, "r") as f:
        expected = f.read()
    assert output == expected


def test_get_molecule(api):
    mol903922 = api.get_molecule(903922)
    assert mol903922.residue_name == "6Y3V"

    assert len(mol903922.atoms) == 45
    assert len(mol903922.bonds) == 44
    assert len(mol903922.angles) == 80
    assert len(mol903922.dihedrals) == 89
    assert len([x for x in mol903922.dihedrals if x.essential]) == 17
    assert len(mol903922.rings) == 0


def test_to_rdkit(mol903922):
    _ = pytest.importorskip("rdkit")
    from rdkit import Chem

    rdmol = mol903922.to_rdkit()
    assert rdmol.GetNumAtoms() == 45
    assert rdmol.GetNumBonds() == 44
    assert Chem.MolToSmarts(rdmol) == (
        "[#6:38](-[H:39])(-[H:40])(-[H:41])-[#6:36](-[H:37])"
        "(-[#6:42](-[H:43])(-[H:44])-[H:45])-[#6:34](=[#8:35])"
        "-[#8:33]-[#6:30](-[H:31])(-[H:32])-[#6:27](-[H:28])(-[H:29])"
        "-[#8:26]-[#6:23](-[H:24])(-[H:25])-[#6:20](-[H:21])(-[H:22])"
        "-[#8:19]-[#6:16](-[H:17])(-[H:18])-[#6:13](-[H:14])(-[H:15])"
        "-[#8:12]-[#6:9](-[H:10])(-[H:11])-[#6:6](-[H:7])(-[H:8])"
        "-[#8:5]-[#6:2](-[H:3])(-[H:4])-[H:1]"
    )


def test_to_mdanalysis(mol903922):
    mda = pytest.importorskip("MDAnalysis")
    u = mol903922.to_mdanalysis()
    assert len(u.atoms) == 45
    assert len(u.bonds) == 44
    assert len(u.angles) == 80
    assert len(u.dihedrals) == 89


@pytest.mark.parametrize("use_input_order, reference", [
    (False, MOL_903922_ITP_AA),
    (True, MOL_903922_ITP_AA_INPUT_ORDER),
])
def test_to_itp(atbmol903922, use_input_order, reference):
    itp_str = atbmol903922.to_itp_string(use_input_order=use_input_order, united=False)
    with open(reference, "r") as f:
        expected = f.read()
    n_skip_lines = 6
    expected = "\n".join(expected.split("\n")[n_skip_lines:])
    itp_str = "\n".join(itp_str.split("\n")[n_skip_lines:])
    assert itp_str == expected


@pytest.mark.parametrize("use_input_order, reference", [
    (False, MOL_903922_ITP_UA),
    (True, MOL_903922_ITP_UA_INPUT_ORDER)
])
def test_to_itp_united(atbmol903922, use_input_order, reference):
    itp_str = atbmol903922.to_itp_string(use_input_order=use_input_order, united=True)
    with open(reference, "r") as f:
        expected = f.read()
    n_skip_lines = 6
    expected = "\n".join(expected.split("\n")[n_skip_lines:])
    itp_str = "\n".join(itp_str.split("\n")[n_skip_lines:])
    assert itp_str == expected
