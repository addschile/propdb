import numpy as np
from typing import List, Optional
from .compute_path import mintree
from .krr import make_krr_weights
from .prop_diabatization import compute_cmat,compute_cmat_krr

def diabatize(nel: int,
              ngeoms: int,
              geoms: np.ndarray,
              eads: np.ndarray,
              Fs: np.ndarray,
              krrflag: Optional[bool] = False,
              krrargs: Optional[List] = None,
              outargs: Optional[List[str]] = None):
  """Computes diabatic energies and transformation matrices from adiabatic
  energies and derivative couplings using propagation diabatization.
  """

  # do some checks about the order of things

  # compute nearest-neighbor path
  nns = mintree(ngeoms, geoms)

  # initialize adiabatic-to-diabatic rotation matrices
  Cmats = np.zeros((ngeoms,nel,nel))
  Cmats[0] = np.eye(nel)

  if krrflag:
    #raise NotImplementedError
    # do diabatization with krr
    nmodes = krrargs[0]
    modes  = krrargs[1]
    if len(krrargs) == 2:
      # default number of integration steps
      nsteps    = 40
      # default krr parameters
      alpha     = 0.5
      gam       = 1.e-4
      modetypes = ['gaussian' for i in range(nmodes)]
    elif len(krrargs) == 6:
      nsteps    = krrargs[2]
      alpha     = krrargs[3]
      gam       = krrargs[4]
      modetypes = krrargs[5]
    else:
      raise ValueError("Incorrect number of arguments for krrargs--must be 2 or 6 arguments.")

    # find krr weights for adiabatic energies
    krr_weights = make_krr_weights(nel,nmodes,ngeoms,modes,eads,alpha,gam,modetypes)

    # do diabatization
    for i in range(1,ngeoms):
      compute_cmat_krr(nns[i], i, nel, nmodes, nns, Cmats, geoms, eads, Fs,
                       nsteps, modes, modetypes, alpha, krr_weights)

  else:

    # do diabatization without krr
    for i in range(1,ngeoms):
      compute_cmat(nns[i], i, nel, nns, Cmats, geoms, eads, Fs)
  
  # output stuff
  if outargs is None:
    dbfilename = 'db.out'
    transfilename = 'trans.out'
  else:
    assert len(outargs) == 2 , "Number of output args must be 2."
    dbfilename = outargs[0]
    transfilename = outargs[1]

  # write information to files
  f = open(dbfilename,'w')
  fmat = open(transfilename,'w')
  for i in range(ngeoms):
    # transform adiabatic potential to diabatic
    vmat = np.dot(Cmats[i].T, np.dot(np.diag(np.array([eads[i,j] for j in range(nel)])), Cmats[i]))
    # write diabatic energies
    f.write('#Record No. : %d\n'%(i+1))
    f.write('%.8f '%(vmat[0,0]))
    f.write('%.8f '%(vmat[0,1]))
    f.write('%.8f\n'%(vmat[1,1]))
    f.flush()
    # write transformation matrix
    fmat.write('#Record No. : %d\n'%(i+1))
    fmat.write('%.8f %.8f\n'%(Cmats[i,0,0],Cmats[i,0,1]))
    fmat.write('%.8f %.8f\n'%(Cmats[i,1,0],Cmats[i,1,1]))
    fmat.flush()
  f.close()
  fmat.close()
