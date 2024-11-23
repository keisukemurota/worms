"""This module contains the functions to generate the Hamiltonian of the Majundar-Ghosh model.

More generally, this model is part of J1-J2-J3 (Triangular) model. When J1=1, J2=J3 = 1/2, the model is the Majundar-Ghosh model.
sps = 2 ** 2 = 4
"""
import numpy as np
from ..core.paulis.spin import Sx, Sz, I2, SxSx, SySy, SzSz
from .. import utils
import logging
from typing import List, Any, Tuple, Dict
from numpy._typing import NDArray
unitary_algorithms = ["original", "2site", "3siteDiag"]


def local(params: Dict[str, Any], D: int = 1) -> Tuple[List[NDArray[Any]], int]:
    """Generate the local Hamiltonian of the Majundar-Ghosh model.

    Available lt is 2, and in this case, the system become a 1D chain.

    Originally, MG model is spin 1/2 model, but here we only consider lt > 2, which means the system is larger than spin 1.
    sps > 4
    """
    if D != 1:
        raise ValueError("J1-J2-J3 is 1D model")
    J1 = params["J1"]
    J2 = params["J2"]
    J3 = params["J3"]
    lt = params["lt"]
    h_bond = SzSz + SxSx + SySy

    u_dimer = np.array([
        [1 / np.sqrt(2), 0, 0, 1 / np.sqrt(2)],
        [1 / np.sqrt(2), 0, 0, -1 / np.sqrt(2)],
        [0, 1 / np.sqrt(2), -1 / np.sqrt(2), 0],
        [0, 1 / np.sqrt(2), 1 / np.sqrt(2), 0],
    ])
    Ud = np.kron(u_dimer, u_dimer)
    if lt == 2:
        # n: MG model has 3 bonds per site.
        _h = utils.sum_ham(h_bond, [[0, 2]], 4, 2) * J1
        _h += utils.sum_ham(h_bond, [[1, 3]], 4, 2) * J1
        _h += utils.sum_ham(h_bond, [[1, 2]], 4, 2) * J3
        _h += utils.sum_ham(h_bond, [[0, 1], [2, 3]], 4, 2) * J2 / 2
        h = _h
        sps = 4
    elif lt == -2:
        # n: MG model has 3 bonds per site.
        _h = utils.sum_ham(h_bond, [[0, 2]], 4, 2) * J1
        _h += utils.sum_ham(h_bond, [[1, 3]], 4, 2) * J1
        _h += utils.sum_ham(h_bond, [[1, 2]], 4, 2) * J3
        _h += utils.sum_ham(h_bond, [[0, 1], [2, 3]], 4, 2) * J2 / 2
        h = Ud @ _h @ Ud.T
        sps = 4
    else:
        raise ValueError("lt != 2 is not implemented")

    return [h], sps


def system(_L: list[int], params: dict) -> Tuple[NDArray[Any], int]:
    """Generate the Hamiltonian of the Majundar-Ghosh model."""
    if len(_L) == 1:
        L = _L[0]
        P = {
            "J1": params["J1"],
            "J2": params["J2"],
            "J3": params["J3"],
            "lt": params["lt"],
        }
        logging.info(f"params : {P}")
        logging.info(f"L      : {L}")
        H_list, sps = local(params, D=1)
        if params["obc"]:
            bonds = [[i, i + 1] for i in range(L - 1)]
        else:
            bonds = [[i, (i + 1) % L] for i in range(L)]

        if len(H_list) != 1:
            raise RuntimeError("something wrong")

        _H = utils.sum_ham(H_list[0], bonds, L, sps)
        return _H, sps

    else:
        raise ValueError("J1-J2-J3(MG) is 1D model")
