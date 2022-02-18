from typing import Union, Dict, Any, Optional
from typing_extensions import Literal

import os
import requests
import pathlib

from .exceptions import InvalidTokenError, ATBQuotaError
from .molecule import Molecule

MOLECULE_TYPES = [
    "heteromolecule",
    "amino acid",
    "nucleic acid",
    "sugar",
    "lipid",
    "solvent",
]


class ATBApi:
    """ATB API

    Examples
    --------

    Create an instance by passing an API token, or the filename of one::

        api = ATBApi(api_token="MY_TOKEN")

    Submit a molecule by passing in a PDB string or filename::

        molid = api.submit_molecule(my_pdb_file.pdb, net_charge=0, molecule_type="heteromolecule")
        assert isinstance(molid, int)

    Get a molecule with a molecule ID::

        molecule = api.get_molecule(molid=molid)
        print(molecule.atoms)
        print(molecule.bonds[0].code)

    Download a molecule file with a molecule ID::

        pdb_as_str = api.download_molecule(molid=903922, format="pdb",
                                           resolution="all_atom", optimized=True)


    """

    host: str = "https://atb.uq.edu.au"
    timeout: int = 45

    def __init__(
        self,
        api_token: Union[str, pathlib.Path] = "~/.ATB_API_TOKEN",
    ):
        api_token_file = pathlib.Path(api_token).expanduser()
        if api_token_file.exists():
            with api_token_file.open("r") as f:
                api_token = f.read()

        self.api_token = api_token
        self.validate_token()

    def validate_token(self):
        self.get_molecule(5441)

    def _check_response(self, response):
        if response.ok:
            return

        if "Could not authenticate" in response.text:
            raise InvalidTokenError(
                f"Invalid authentication credentials. Using api_token={self.api_token}"
            )
        contents = self.decode_yaml(response.text)
        if "daily quota" in contents.get("error", ""):
            raise ATBQuotaError(contents["error"])

    def get_url(self, namespace: str, endpoint: str):
        if not endpoint.endswith(".py"):
            endpoint = endpoint + ".py"
        return os.path.join(self.host, "api", "current", namespace, endpoint)

    @staticmethod
    def decode_yaml(text: str) -> Dict[Any, Any]:
        import yaml

        return yaml.load(text, yaml.Loader)

    def get(
        self,
        namespace: str,
        endpoint: str,
        data: Dict[str, Any] = {},
        decode: bool = True,
    ) -> Union[str, Dict[Any, Any]]:
        url = self.get_url(namespace, endpoint)
        params = {
            "api_token": self.api_token,
            "api_format": "yaml",
        }
        params.update(data)
        response = requests.get(url, params=params)
        self._check_response(response)
        text = response.text
        if decode:
            text = self.decode_yaml(text)
        return text

    def post(
        self,
        namespace: str,
        endpoint: str,
        data: Dict[str, Any] = {},
    ):
        url = self.get_url(namespace, endpoint)
        params = {"api_token": self.api_token}
        params.update(data)
        response = requests.post(url, data=params)
        self._check_response(response)
        return self.decode_yaml(response.text)

    def download_molecule_file(
        self,
        molid: int,
        format="itp",
        resolution="united_atom",
        optimized: bool = True,
        filename: Optional[str] = None,
    ):
        FORMAT = {
            "pdb": "pdb",
            "mtb": "mtb",
            "itp": "rtp",
        }
        RESOLUTION = {
            "united_atom": "uniatom",
            "ua": "uniatom",
            "all_atom": "allatom",
            "aa": "allatom",
        }
        fmt = FORMAT.get(format.lower(), format)
        res = RESOLUTION.get(resolution.lower(), resolution)

        if fmt == "pdb":
            if optimized:
                res = f"{res}_optimised"
            else:
                res = f"{res}_unoptimised"

        mol_params = dict(
            molid=molid,
            hash="LATEST",
            outputType="top",
            ffVersion="54A7",
            file=f"{fmt}_{res}",
        )

        text = self.get(
            "molecules",
            "download_file.py",
            data=mol_params,
            decode=False,
        )
        if filename is not None:
            with open(str(filename), "w") as f:
                f.write(text)
        return text

    def submit_molecule(
        self,
        pdb: str,
        net_charge: int = 0,
        molecule_type: Literal[(*MOLECULE_TYPES,)] = "heteromolecule",  # type: ignore
        public: bool = True,
    ) -> int:
        from .utils import read_from_string_or_file

        pdb = read_from_string_or_file(pdb)

        data = dict(
            pdb=pdb,
            netcharge=net_charge,
            moltype=molecule_type,
            public=public,
        )
        text = self.post("molecules", "submit.py", data)
        if "error_msg" in text:
            molid = int(text["error_msg"].split("molid=")[-1].split(",")[0])
        else:
            molid = text["molid"]
        return molid

    def get_molecule(self, molid: int):
        info = self.get(
            "molecules",
            "generate_mol_data",
            {"molid": molid},
        )
        return Molecule.from_atb_dict(info)
