import numpy as np
from typing import List , Optional

def make_gweights(nmodes: int,
                  alpha: float,
                  coord: np.ndarray,
                  modes: np.ndarray,
                  modetypes: List[str]) -> np.ndarray:
  """Evaluate the gaussian weight vector at a single point.
  """
  gweights = np.zeros(modes.shape[1])
  for i in range(nmodes):
    if modetypes[i] == 'gaussian': 
      gweights += (coord[i]-modes[i])**2.
    elif modetypes[i] == 'periodic': 
      gweights += 2.*(np.sin(0.5*(coord[i] - modes[i])))**2.
  return np.exp(-alpha*gweights)

def make_gweights_mat(nmodes: int,
                      ngeoms: int,
                      alpha: float,
                      modes: np.ndarray,
                      modetypes: List[str],
                      gam: Optional[float] = 0.0) -> np.ndarray:
  """Evaluate the gaussian weight matrix at all points.
  """
  gweights = np.zeros((ngeoms,ngeoms))
  for i in range(nmodes):
    if modetypes[i] == 'gaussian': 
      gweights += (np.outer(modes[i],np.ones(ngeoms)) 
                      - np.outer(np.ones(ngeoms),modes[i]))**2
    elif modetypes[i] == 'periodic': 
      gweights += 2.*(np.sin(0.5*(np.outer(modes[i],np.ones(ngeoms)) 
                      - np.outer(np.ones(ngeoms),modes[i]))))**2
  return np.exp(-alpha*gweights + gam**2.*np.eye(gweights.shape[0]))

def make_vkrr(nmodes: int,
              ngeoms: int,
              alpha: float,
              modes: np.ndarray,
              modetypes: List[str],
              krr_weights: np.ndarray) -> np.ndarray:
  """Computes the energy at all points from kernel ridge regression potential 
  by contrating the gaussian kernel with the krr weights.
  """
  gmat = make_gweights_mat(nmodes,ngeoms,alpha,modes,modetypes)
  return np.tensordot(krr_weights,gmat,axes=[[0],[1]])

def krr_energy(nmodes: int,
               alpha: float,
               coord: np.ndarray,
               modes: np.ndarray,
               modetypes: List[str],
               krr_weights: np.ndarray) -> float:
  """Computes the energy at a single point from kernel ridge regression 
  potential by contrating the gaussian kernel with the krr weights.
  """
  gs = make_gweights(nmodes,alpha,coord,modes,modetypes)
  return np.tensordot(krr_weights,gs,axes=[[0],[0]])

def make_krr_weights(nel: int,
                     nmodes: int,
                     ngeoms: int,
                     modes: np.ndarray,
                     eads: np.ndarray,
                     alpha: float,
                     gam: float,
                     modetypes: List[str]) -> np.ndarray:
  """Computes the weights obtained from kernel-ridge regression."""
  Kmat = make_gweights_mat(nmodes,ngeoms,alpha,modes,modetypes,gam=gam)
  Kinv = np.linalg.inv(Kmat)
  krr_weights = np.zeros((ngeoms,nel))
  for i in range(nel):
    krr_weights[:,i] = np.dot(Kinv,eads[:,i])
  return krr_weights
