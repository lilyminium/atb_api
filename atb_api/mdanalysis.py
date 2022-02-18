from .utils import require_package  # noqa

require_package("MDAnalysis")  # noqa

from MDAnalysis.core.topologyattrs import AtomAttr
import numpy as np


class PartialCharge(AtomAttr):
    attrname = "partial_charges"
    singular = "partial_charge"

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na)


class United(AtomAttr):
    attrname = "united"
    singular = "united"

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.zeros(na, dtype=bool)


class OutputAtomisticID(AtomAttr):
    attrname = "output_atomistic_ids"
    singular = "output_atomistic_id"
    dtype = int

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.arange(1, na + 1)


class OutputUnitedID(AtomAttr):
    attrname = "output_united_ids"
    singular = "output_united_id"
    dtype = int

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.arange(1, na + 1)
