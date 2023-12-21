import torch
import numpy as np
from torch import nn

# import math
from typing import Union, List
import logging


class MinimumEnergyLoss(nn.Module):
    def __init__(
        self,
        h_list: List[Union[np.ndarray, torch.Tensor]],
        device: torch.device = torch.device("cpu"),
        decay: float = 0.1,
        # n : after "decay" step, the regularization term will e^(-1) times smaller
    ):
        super(MinimumEnergyLoss, self).__init__()

        offset = []
        self.shift_origin_offset = []
        self.h_list = []
        self.X = []
        self.weight = 1 if decay > 0 else 0
        if decay < 0:
            raise ValueError("decay should be positive.")
        self.weight_decay = decay
        for i in range(len(h_list)):
            if isinstance(h_list[i], np.ndarray):
                self.h_list.append(torch.from_numpy(
                    h_list[i]).to(torch.float64).to(device))
            elif isinstance(h_list[i], torch.Tensor):
                self.h_list.append(h_list[i].to(device))
            else:
                raise TypeError(
                    "h should be of type np.ndarray or torch.Tensor.")
            E, V = torch.linalg.eigh(self.h_list[i])

            logging.info(f"maximum energy of local hamiltonian {i}: {E[-1]:.3f}")
            logging.info(f"minimum energy of local hamiltonian {i}: {E[0]:.3f}")
            offset.append(E[-1])
            self.h_list[i][:] = self.h_list[i] - offset[i] * \
                torch.eye(h_list[i].shape[1], device=device)
            logging.info(f"offset of local hamiltonian {i}: {offset[i]:.3f}")
            self.X.append(V[:, 0].to(device))
            self.shift_origin_offset.append(offset[i] - E[0])

    def forward(self, U: torch.Tensor) -> torch.Tensor:
        """Return loss value for the minimum eigen loss. (Or minimum eigen local loss).

        Local Hamiltonian will be shifted so that miminum value will be 0
        """
        # add minimum energy of each local hamiltonian for calculating the total energy
        # initialize with torch tensor
        loss = torch.zeros(1, device=U.device)
        for i in range(len(self.h_list)):
            loss += self.minimum_energy_loss(self.h_list[i], U) - self.shift_origin_offset[i]

        self.weight = self.weight * np.exp(-1 / self.weight_decay)
        return torch.abs(loss)

    def minimum_energy_loss(self, H: torch.Tensor, U: torch.Tensor) -> torch.Tensor:
        """Calculate the minimum energy of a system using the reverse iteration method.

        The first ground state is calculated using the eigendecomposition of the Hamiltonian.
        """
        A = U @ H @ U.T
        result_abs = self.get_stoquastic(A)

        # n: corresponds to caluclate eneergy of maximum superposition state
        negativity = torch.abs(A - result_abs).mean() / H.shape[0]
        try:
            E = torch.linalg.eigvalsh(result_abs)
        except RuntimeError:
            # If there are some errors during the eigen decomposition.
            result_abs = (result_abs + result_abs.T)/2
            E = torch.linalg.eigvalsh(result_abs)

        return - E[0] + self.weight * negativity

    def stoquastic(self, A: torch.Tensor):
        """Change the sign of all non-diagonal elements into negative.

        If A is already negative definite matrix,
        then just taking negative absolute value sufficies.
        """
        return -torch.abs(A)

    def get_stoquastic(self, A: torch.Tensor) -> torch.Tensor:
        """Return the stoquastic matrix of a given matrix.

        First, calculate the matrix A = U @ h @ U.T. Then apply stoquastic function.
        """
        return self.stoquastic(A)


class SystemQUasiEnergyLoss(MinimumEnergyLoss):
    def __init__(
        self,
        h_list: List[Union[np.ndarray, torch.Tensor]],
        device: torch.device = torch.device("cpu"),
    ):
        # h_list must be 1 element list
        if len(h_list) != 1:
            raise ValueError("There should be only one local Hamiltonian for SystemQuasiEnergyLoss")

        super(SystemQUasiEnergyLoss, self).__init__(h_list, device)

    def forward(self, U: torch.Tensor, nr: int = 2, regularizer: float = 0.0) -> torch.Tensor:
        """
        Return loss value for the System Quasi Energy Loss.
        Quasi Energy loss is approximated value of System Energy Loss.
        """
        if nr < 1:
            raise ValueError("nr should be positive integer.")

        X = self.X[0]
        H = self.h_list[0]
        offset = self.shift_origin_offset[0]
        return self.system_quasi_energy_loss(U, X, H, nr, regularizer) - offset

    def system_quasi_energy_loss(
            self,
            U: torch.Tensor,
            X: torch.Tensor,
            H: torch.Tensor,
            nr: int = 2,
            regularizer: float = 0.0) -> torch.Tensor:
        """
        Return loss value for the System Quasi Energy Loss.
        Quasi Energy loss is approximated value of System Energy Loss.
        """
        if nr < 1:
            raise ValueError("nr should be positive integer.")

        SUx = torch.abs(U @ X)
        SUH = self.get_stoquastic(U @ H @ U.T)
        y = SUx
        for _ in range(10):
            y = SUH @ y
            if _ % 4 == 0:
                y = y / torch.norm(y)

        quasi_Sgs = torch.abs(y / torch.norm(y))
        z = SUH @ quasi_Sgs
        res = -(quasi_Sgs @ z)
        res += regularizer * (1 - torch.abs(quasi_Sgs.dot(z) / torch.norm(z)))
        return res

        # class SystemEnergyLoss(nn.Module):
        #     def __init__(
        #         self,
        #         H: Union[np.ndarray, torch.Tensor],
        #         device: torch.device = torch.device("cpu"),
        #         P_single: Union[np.ndarray, torch.Tensor] = None,
        #     ):
        #         super(SystemEnergyLoss, self).__init__()
        #
        #         if isinstance(H, np.ndarray):
        #             self.H = torch.from_numpy(H).to(device)
        #         elif isinstance(H, torch.Tensor):
        #             self.H = H.to(device)
        #         else:
        #             raise TypeError("H should be of type np.ndarray or torch.Tensor.")
        #
        #         if self.H.ndim == 2:
        #             E, V = torch.linalg.eigh(self.H)
        #         elif self.H.ndim == 3:
        #             # sum over the first axis
        #             E, V = torch.linalg.eigh(self.H.sum(axis=0))
        #             if P_single is None:
        #                 raise ValueError("P should be provided for 3D H.")
        #             if P_single.shape != self.H[0].shape:
        #                 raise ValueError("P should have the same shape as H[0].")
        #             self.P = torch.from_numpy(P_single).to(device)
        #             self.P_prime = torch.from_numpy(1 - P_single).to(device)
        #         offset = E[-1]
        #         self.eye_offset = offset * \
        #             torch.eye(self.H.shape[1], device=device)
        #         self.X = V[:, 0].to(device)
        #
        #     def forward(self, U: torch.Tensor) -> torch.Tensor:
        #         # A = torch.matmul(U, torch.matmul(self.H, U.T))
        #         # result_abs = self.stoquastic(A)
        #         result_abs = self.get_stoquastic(self.H, U)
        #         E = torch.linalg.eigvalsh(result_abs)
        #         z = torch.exp(-E * 1).sum()
        #         return torch.log(z)
        #
        #     def stoquastic(self, A: torch.Tensor):
        #         return -torch.abs(A - self.eye_offset) + self.eye_offset
        #
        #     def initializer(self, U: Union[torch.Tensor, None] = None):
        #         raise NotImplementedError(
        #             "Initializer is not implemented for QuasiEnergyLoss.")
        #         return
        #
        #     def get_stoquastic(self, H: torch.Tensor, U: torch.Tensor) -> torch.Tensor:
        #         """
        #         Return the stoquastic matrix of a given matrix.
        #         """
        #         A = U @ H @ U.T
        #         if A.ndim == 2:
        #             return self.stoquastic(A)
        #         elif A.ndim == 3:
        #             tmp1 = self.P * self.stoquastic((A.sum(axis=0)))
        #             tmp2 = self.P_prime * self.stoquastic(A).sum(axis=0)
        #             return tmp1 + tmp2
        #
        #
        # class SystemQuasiEnergyLoss(SystemEnergyLoss):
        #     def __init__(
        #         self,
        #         H: Union[np.ndarray, torch.Tensor],
        #         N: int = 10,
        #         r: float = 0,  # * regularization
        #         device: torch.device = torch.device("cpu"),
        #         P_single: Union[np.ndarray, torch.Tensor] = None,
        #     ):
        #         super(SystemQuasiEnergyLoss, self).__init__(H, device, P_single)
        #         if self.H.ndim == 2:
        #             self.H = self.H - self.eye_offset
        #         elif self.H.ndim == 3:
        #             self.H = self.H - self.eye_offset / self.H.shape[0]
        #         self.N = int(N)
        #
        #     def forward(self, U: torch.Tensor, r: float = 0) -> torch.Tensor:
        #         SUx = torch.abs(U @ self.X)
        #         SUH = self.get_stoquastic(self.H, U)
        #         y = SUx
        #         for _ in range(self.N):
        #             y = SUH @ y
        #             if _ % 4 == 0:
        #                 y = y / torch.norm(y)
        #         quasi_Sgs = torch.abs(y / torch.norm(y))
        #
        #         z = SUH @ quasi_Sgs
        #         return -(quasi_Sgs @ z + offset) + r * \
        #             (1 - torch.abs(quasi_Sgs.dot(z) / torch.norm(z)))
        #
        #     def initializer(self, U: Union[torch.Tensor, None] = None):
        #         raise NotImplementedError(
        #             "Initializer is not implemented for QuasiEnergyLoss.")
        #         return
        #
        #
        # class SystemMinimumEnergyLoss(nn.Module):
        #     """
        #     Calculate the minimum energy of a system using the reverse iteration method.
        #     The first ground state is calculated using the eigendecomposition of the Hamiltonian.
        #     The next ground states will be calculated using the reverse iteration method.
        #     You need to use small learning rates for this loss function in order to converge.
        #     """
        #
        #     def __init__(
        #         self,
        #         H: Union[np.ndarray, torch.Tensor],
        #         device: torch.device = torch.device("cpu"),
        #         P_single: Union[np.ndarray, torch.Tensor] = None,
        #     ):
        #         super(SystemMinimumEnergyLoss, self).__init__()
        #
        #         if isinstance(H, np.ndarray):
        #             self.H = torch.from_numpy(H).to(device)
        #         elif isinstance(H, torch.Tensor):
        #             self.H = H.to(device)
        #         else:
        #             raise TypeError("H should be of type np.ndarray or torch.Tensor.")
        #         if self.H.ndim == 2:
        #             E, V = torch.linalg.eigh(self.H)
        #         elif self.H.ndim == 3:
        #             # sum over the first axis
        #             E, V = torch.linalg.eigh(self.H.sum(axis=0))
        #             if P_single is None:
        #                 raise ValueError("P should be provided for 3D H.")
        #             if P_single.shape != self.H[0].shape:
        #                 raise ValueError("P should have the same shape as H[0].")
        #             self.P = torch.from_numpy(P_single).to(device)
        #             self.P_prime = torch.from_numpy(1 - P_single).to(device)
        #
        #         offset = E[-1]
        #         self.eye = torch.eye(self.H.shape[1], device=device)
        #         self.eye_offset = offset * self.eye
        #         self.X = V[:, 0].to(device)
        #         self.V_old = None  # * V_tilde
        #         self.E_old = None  # * E_min_tilde
        #
        #     def forward(self, U: torch.Tensor) -> torch.Tensor:
        #         H_tilde = self.get_stoquastic(self.H, U)
        #         with torch.no_grad():
        #             self.V_old = torch.linalg.solve(
        #                 H_tilde - self.eye * self.E_old, self.V_old)
        #             self.V_old = self.V_old / torch.norm(self.V_old)
        #         E = self.V_old.detach() @ H_tilde @ self.V_old.detach()
        #         self.E_old = E
        #         return -E
        #
        #     def stoquastic(self, A: torch.Tensor):
        #         return -torch.abs(A - self.eye_offset) + self.eye_offset
        #
        #     def initializer(self, U: Union[torch.Tensor, None] = None):
        #         if U is None:
        #             U = torch.eye(
        #                 self.H.shape[1], device=self.H.device, dtype=self.H.dtype)
        #         U = U.detach()
        #         H_tilde = self.get_stoquastic(self.H, U)
        #         E, V = torch.linalg.eigh(H_tilde)
        #         self.V_old = V[:, 0]
        #         self.E_old = E[0]
        #         return
        #
        #     def get_stoquastic(self, H: torch.Tensor, U: torch.Tensor) -> torch.Tensor:
        #         """
        #         Return the stoquastic matrix of a given matrix.
        #         """
        #         A = U @ H @ U.T
        #         if A.ndim == 2:
        #             return self.stoquastic(A)
        #         elif A.ndim == 3:
        #             tmp1 = self.P * self.stoquastic((A.sum(axis=0)))
        #             tmp2 = self.P_prime * self.stoquastic(A).sum(axis=0)
        #             return tmp1 + tmp2
        #             # return self.stoquastic(A).sum(axis=0)
