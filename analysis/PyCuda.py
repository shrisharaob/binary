import pycuda.gpuarray as gpuarray
import pycuda.autoinit
import numpy as np
from skcuda import linalg
#linalg.init()
# Compute right eigenvectors of a symmetric matrix A and verify A*vr = vr*w
a = np.array(([1,3],[3,5]), np.float32, order='F')
a_gpu = gpuarray.to_gpu(a)
vr_gpu, w_gpu = linalg.eig(a_gpu, 'N', 'V')
np.allclose(np.dot(a, vr_gpu.get()), np.dot(vr_gpu.get(), np.diag(w_gpu.get())), 1e-4)
# Compute left eigenvectors of a symmetric matrix A and verify vl.T*A = w*vl.T
a = np.array(([1,3],[3,5]), np.float32, order='F')
a_gpu = gpuarray.to_gpu(a)
w_gpu, vl_gpu = linalg.eig(a_gpu, 'V', 'N')
np.allclose(np.dot(vl_gpu.get().T, a), np.dot(np.diag(w_gpu.get()), vl_gpu.get().T), 1e-4)
# Compute left/right eigenvectors of a symmetric matrix A and verify A = vr*w*vl.T
a = np.array(([1,3],[3,5]), np.float32, order='F')
a_gpu = gpuarray.to_gpu(a)
vr_gpu, w_gpu, vl_gpu = linalg.eig(a_gpu, 'V', 'V')
np.allclose(a, np.dot(vr_gpu.get(), np.dot(np.diag(w_gpu.get()), vl_gpu.get().T)), 1e-4)
# Compute eigenvalues of a square matrix A and verify that trace(A)=sum(w)
a = np.array(np.random.rand(9,9), np.float32, order='F')
a_gpu = gpuarray.to_gpu(a)
w_gpu = linalg.eig(a_gpu, 'N', 'N')
np.allclose(np.trace(a), sum(w_gpu.get()), 1e-4)
# Compute eigenvalues of a real valued matrix A possessing complex e-valuesand
a = np.array(np.array(([1, -2], [1, 3])), np.float32, order='F')
a_gpu = gpuarray.to_gpu(a)
w_gpu = linalg.eig(a_gpu, 'N', 'N', imag='T')
# Compute eigenvalues of a complex valued matrix A and verify that trace(A)=sum(w)
a = np.array(np.random.rand(2,2) + 1j*np.random.rand(2,2), np.complex64, order='F')
a_gpu = gpuarray.to_gpu(a)
w_gpu = linalg.eig(a_gpu, 'N', 'N')
np.allclose(np.trace(a), sum(w_gpu.get()), 1e-4)
