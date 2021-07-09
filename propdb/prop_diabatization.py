import numpy as np
from scipy.linalg import expm
from typing import List
from .krr import krr_energy

def update_Ftrans(nel,eads,Fs,dR,Ftrans):
  """
  """
  count = 0
  for i in range(nel-1):
    for j in range(i+1,nel):
      Ftrans[i,j] += 0.5 * np.sum(Fs[count]*dR) / (eads[j]-eads[i])
      count += 1

def antihermitize(nel,F):
  """
  """
  for i in range(nel-1):
    for j in range(i+1,nel):
      F[j,i] = -F[i,j]

def compute_cmat(ind0: int,
                 indf: int,
                 nel: int,
                 nns: List[int],
                 Cmats: np.ndarray,
                 geoms: np.ndarray,
                 eads: np.ndarray,
                 Fs: np.ndarray):
  """Compute adiabatic-to-diabatic rotation matrix using propagation between two
  geometries. If rotation matrix for starting geometry hasn't been computed,
  this function gets recursively called until reaching a geometry for which the
  rotation matrix has been computed. Calculates line integrals using trapezoid
  rule approximation of integral over derivative coupling between two geometries.

  Parameters
  ----------
  ind0: int, index of starting geometry in line integral
  indf: int, index of ending geometry in line integral
  nel: int, number of electronic states
  nns: list, indexes of nearest neighbors starting from an origin
  Cmats: list, adiabatic-to-diabatic transformation matrices
  geoms: np.ndarray, numpy array containing geometries as numpy arrays
  eads: np.ndarray, numpy array containing adiabatic energies at each geometry
  Fs: np.ndarray, numpy array containing nonadiabatic coupling vectors as numpy arrays
  """

  if not Cmats[indf].any():
    # compute current Cmat
    if not Cmats[ind0].any():
      # compute previous Cmat
      compute_cmat(nns[ind0],ind0,nel,nns,Cmats,geoms,eads,Fs)

    # compute current Cmat from previous Cmat
    Fprev = Fs[ind0].copy()
    Fnext = Fs[indf].copy()

    # compute phase information
    for i in range(int((nel*nel-nel)/2)):
      sgn = np.sum((Fprev[i]*Fnext[i]))
      if sgn < 0.0:
        Fnext[i] *= -1.0

    # prep line integral
    Ftrans = np.zeros((nel,nel))

    # trapezoid rule integration
    count = 0
    for i in range(nel-1):
      for j in range(i+1,nel):
        Ftrans[i,j] += 0.5*np.sum((Fprev[count] + Fnext[count])*(geoms[indf]-geoms[ind0]))
        Ftrans[j,i] = -Ftrans[i,j]
        count += 1

    # make symmetric propagator
    try:
      CpropR = expm(-0.5*Ftrans)
      CpropL = expm(0.5*Ftrans)
      Cprop = np.dot(np.linalg.inv(CpropL),CpropR)
    except:
      print('symmetric propagator construction failed')

    # propagate
    Cmats[indf] = np.dot(Cprop,Cmats[ind0])

    return

  else:

    return

def compute_cmat_krr(ind0: int,
                     indf: int,
                     nel: int,
                     nmodes: int,
                     nns: List[int],
                     Cmats: np.ndarray,
                     geoms: np.ndarray,
                     eads: np.ndarray,
                     Fs: np.ndarray,
                     nsteps: int,
                     modes: np.ndarray,
                     modetypes: List[str],
                     alpha: float,
                     krr_weights: np.ndarray):
  """Compute adiabatic-to-diabatic rotation matrix using propagation between two
  geometries. If rotation matrix for starting geometry hasn't been computed,
  this function gets recursively called until reaching a geometry for which the
  rotation matrix has been computed. Calculates line integrals using trapezoid
  rule approximation of integral over derivative coupling divided by adiabatic 
  energy between two geometries. Adiabatic energies are approximated by kernel
  ridge regression previously calculated.

  Parameters
  ----------
  ind0: int, index of starting geometry in line integral
  indf: int, index of ending geometry in line integral
  nel: int, number of electronic states
  nmodes: int, number of modes for krr 
  nns: list, indexes of nearest neighbors starting from an origin
  Cmats: list, adiabatic-to-diabatic transformation matrices
  geoms: np.ndarray, numpy array containing geometries as numpy arrays
  eads: np.ndarray, numpy array containing adiabatic energies at each geometry
  Fs: np.ndarray, numpy array containing nonadiabatic coupling vectors as numpy arrays
  nsteps: int, number of steps used for trapezoid rule integration
  modes: np.ndarray, numpy array containing values of the modes
  modetypes: list, types of distance measures for kernel
  alpha: float, parameter giving decay of kernel distance measure
  krr_weights: np.ndarray, kernel weights
  """

  if not Cmats[indf].any():
    # compute current Cmat
    if not Cmats[ind0].any():
      # compute previous Cmat
      compute_cmat_krr(nns[ind0],ind0,nel,nmodes,nns,Cmats,geoms,eads,Fs,
                       nsteps,modes,modetypes,alpha,krr_weights)

    # compute current Cmat from previous Cmat
    Fprev = Fs[ind0].copy()
    Fnext = Fs[indf].copy()

    # compute phase information
    for i in range(int((nel*nel-nel)/2)):
      sgn = np.sum((Fprev[i]*Fnext[i]))
      if sgn < 0.0:
        Fnext[i] *= -1.0

    # prep line integral
    Ftrans = np.zeros((nel,nel))

    # starting positions of integral
    Rnm = modes[:,ind0].copy()

    # differentials along integral
    dR   = (geoms[indf]-geoms[ind0])/float(nsteps)
    dRnm = (modes[:,indf]-modes[:,ind0])/float(nsteps)
    dF   = (Fnext - Fprev)/float(nsteps)

    # trapezoid rule integration
    for j in range(nsteps):
      if j==0:
        eapprox = eads[ind0]
      elif j==nsteps-1:
        eapprox = eads[indf]
      else:
        # compute adiabatic energies from krr
        eapprox = 0.5*krr_energy(nmodes,alpha,Rnm,modes,modetypes,krr_weights)
      update_Ftrans(nel,eapprox,Fprev,dR,Ftrans)
      Rnm += dRnm
      Fprev += dF
    antihermitize(nel,Ftrans)

    try:
      # make symmetric propagator
      CpropR = expm(-0.5*Ftrans)
      CpropL = expm(0.5*Ftrans)
      Cprop = np.dot(np.linalg.inv(CpropL),CpropR)
    except:
      print(Ftrans,CpropR,CpropL)

    # propagate
    Cmats[indf] = np.dot(Cprop,Cmats[ind0])

    return

  else:

    return
