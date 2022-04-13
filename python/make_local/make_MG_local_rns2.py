"""
reduce negative sign with unitary transformation
"""
import numpy as np
from scipy import sparse
import os
import sys
import torch
import torch.optim

sys.path.insert(0, "../") 
from nsp.utils.lossfunc import make_positive_np
from nsp.utils import optm
from functions import *



I = np.identity(2)

Sz = np.zeros([2,2])
Sz[0,0] = 1/2
Sz[1,1] = -1/2
Sx = np.zeros([2,2])
Sx[1,0] = 1/2
Sx[0,1] = 1/2
Sy = np.zeros([2,2], dtype=np.complex128)
Sy[1,0] = 1j/2
Sy[0,1] = -1j/2


SzSz = np.kron(Sz,Sz).astype(np.float64)
SxSx = np.kron(Sx,Sx).astype(np.float64)
SySy = np.kron(Sy,Sy).astype(np.float64)

lh = SzSz + SxSx + SySy
# lh = - lh


LH_ = sparse.csr_matrix((2**3,2**3), dtype = np.float64)
i = 0
LH_ += l2nl(lh/2, 3, [0, 1], sps = 2)
LH_ += l2nl(lh/2, 3, [0, 2], sps = 2)
LH_ += l2nl(lh/2, 3, [1, 2], sps = 2)


LH = sparse.csr_matrix((2**6,2**6), dtype = np.float64)
LH += l2nl(LH_/2, 6, [0, 1, 2], sps = 2)
LH += l2nl(LH_, 6, [1, 2, 3], sps = 2)
LH += l2nl(LH_, 6, [2, 3, 4], sps = 2)
LH += l2nl(LH_/2, 6, [3, 4, 5], sps = 2)


H = sparse.csr_matrix((2**12,2**12), dtype = np.float64)
H += l2nl(LH, 4, [0, 1], sps = 8)
H += l2nl(LH, 4, [1, 2], sps = 8)

# LH2 = l2nl(LH, 2, [1, 0], sps = 8)
X = -H.toarray()
# X1 = -LH.toarray()
# X2 = -LH2.toarray()

Y = make_positive_np(X)


E = np.linalg.eigvalsh(X)[-1]
E_init = np.linalg.eigvalsh(Y)[-1]
print("maximum eigenvalue (target loss) : ", E)
print("initial loss                     : ", E_init)


"""
Dual Annealing
"""

targets_da = []

# def callbackF(x, f, context):
#     if context == 1:
#         print("target value : {:.5f} in the context {}".format(f,context))
#     targets_da.append(f)

# func = optm.unitary_optm([X1, X2], 28, add = True)
# import scipy.optimize as optimize
# bounds = [[-100, 100] for _ in range(28)]
# ret = optimize.dual_annealing(
#     func, bounds = bounds, restart_temp_ratio = 1e-5, visit = 2.7, initial_temp = 3*10**4, maxiter = 500, callback = callbackF)

# func(ret.x)
# X_da = np.zeros_like(X1)
# for mat in [X1, X2]:
#     X_da += make_positive_np(func.U @ mat @ func.U.T)
# E_da = np.linalg.eigvalsh(X_da)[-1]
# print("dual annealing loss               : ", E_da)


# model, gl = optm.optim_matrix_symm(
#     [torch.tensor(X1), torch.tensor(X2)], 100000,
#     optm_method = torch.optim.SGD, seed = 10, 
#     lr = 0.001, add = True)

# U = model.matrix
# X_loss1 = np.zeros_like(X1)
# for mat in [X1, X2]:
#     X_loss1 += make_positive_np(U @ mat @ U.T)
# E_loss1 = np.linalg.eigvalsh(X_loss1)[-1]
# print("dual annealing loss               : ", E_loss1)


import torch.optim
model, gl = optm.optim_matrix_symm(
        [torch.tensor(X)],
        10000, 
        optm_method = optm.scheme1, 
        gamma = 0.001,
        r = 1,
        )
U = model.matrix

X = np.array(model(torch.tensor(X)).data[0], dtype=np.float64)
path = "../array/MG_union_rns2_bond"
if not os.path.isfile(path):
  np.save(path,X)
  print("save : ", path+".npy")
  beauty_array(X,path + ".txt")