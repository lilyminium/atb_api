atb_api
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/lilyminium/atb_api/workflows/CI/badge.svg)](https://github.com/lilyminium/atb_api/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/lilyminium/atb_api/branch/master/graph/badge.svg)](https://codecov.io/gh/lilyminium/atb_api/branch/master)



This API is a re-write of the [UQ version](https://github.com/ATB-UQ/atb_api_public) in Python 2.7

This is a work in progress. It is for Lily's personal use, so not all features of the old API are supported yet.
Please open an issue on the [Issue tracker](https://github.com/lilyminium/atb_api/issues) if you notice a bug or need a feature.


### Example

You will need a valid ATB token to use this API. Please email the ATB administrators to request a token.

Create an instance by passing an API token, or the filename of one:

```python
    api = ATBApi(api_token="MY_TOKEN")
```

Submit a molecule by passing in a PDB string or filename:

```python
    molid = api.submit_molecule(my_pdb_file.pdb, net_charge=0, molecule_type="heteromolecule")
    assert isinstance(molid, int)
```

Get a molecule with a molecule ID:

```python
    molecule = api.get_molecule(molid=molid)
    print(molecule.atoms)
    print(molecule.bonds[0].code)
```

Download a molecule file with a molecule ID:

```python
    pdb_as_str = api.download_molecule(molid=903922, format="pdb",
                                       resolution="all_atom", optimized=True)
```

### Copyright

Copyright (c) 2022, Lily Wang


#### Acknowledgements
 
Many thanks to the original [ATB API](https://github.com/ATB-UQ/atb_api_public) from Bertrand Caron.

Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
