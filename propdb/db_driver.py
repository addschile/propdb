import numpy as np
from time import time
from typing import List, Optional
from .compute_path import mintree
from .krr import make_krr_weights,make_vkrr
from .prop_diabatization import compute_cmat,compute_cmat_krr

def diabatize(nel: int,
              ngeoms: int,
              nmodes: int,
              geoms: np.ndarray,
              modes: np.ndarray,
              eads: np.ndarray,
              Fs: np.ndarray,
              krrflag: Optional[bool] = False,
              krrargs: Optional[List] = None,
              outargs: Optional[List[str]] = None):
  """Computes diabatic energies and transformation matrices from adiabatic
  energies and derivative couplings using propagation diabatization.
  """

  # TODO do some checks about the order of things

  # compute nearest-neighbor path
  nns = mintree(ngeoms, modes.T)

  # initialize adiabatic-to-diabatic rotation matrices
  Cmats = np.zeros((ngeoms,nel,nel))
  Cmats[0] = np.eye(nel)

  if ngeoms>10:
    every = int(ngeoms/10)
  else:
    every = 1

  if krrflag:
    if krrargs is None:
      # default number of integration steps
      nsteps    = 40
      # default krr parameters
      alpha     = 0.5
      gam       = 1.e-4
      modetypes = ['gaussian' for i in range(nmodes)]
    elif len(krrargs) == 5:
      nsteps    = krrargs[0]
      alpha     = krrargs[1]
      gam       = krrargs[2]
      modetypes = krrargs[3]
      vkrrflag  = krrargs[4]
    else:
      raise ValueError("Incorrect number of arguments for krrargs--must be 2 or 6 arguments.")

    # find krr weights for adiabatic energies
    krr_weights = make_krr_weights(nel,nmodes,ngeoms,modes,eads,alpha,gam,modetypes)

    # print the potential from krr
    if vkrrflag:
      vkrr = make_vkrr(nmodes,ngeoms,modes,alpha,modetypes,krr_weights)
      np.save('vkrr.npy',vkrr,allow_pickle=True)

    # do diabatization
    btime = time()
    for i in range(1,ngeoms):
      if i%every==0:
        rtime = time()-btime
        print('%d Percent done.........%.3f'%((i/every)*10,rtime))
      compute_cmat_krr(nns[i], i, nel, nmodes, nns, Cmats, geoms, eads, Fs,
                       nsteps, modes, modetypes, alpha, krr_weights)

  else:

    # do diabatization without krr
    btime = time()
    for i in range(1,ngeoms):
      if i%every==0:
        rtime = time()-btime
        print('%d Percent done.........%.3f'%((i/every)*10,rtime))
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
    f.write('#Geometry %d\n'%(i+1))
    for j in range(nel):
      for k in range(j,nel):
        f.write('%.8f '%(vmat[j,k]))
    f.write('\n')
    f.flush()
    # write transformation matrix
    fmat.write('#Geometry %d\n'%(i+1))
    for j in range(nel):
      for k in range(nel):
        fmat.write('%.8f '%(Cmats[i,j,k]))
      fmat.write('\n')
    fmat.write('\n')
    fmat.flush()
  f.close()
  fmat.close()
