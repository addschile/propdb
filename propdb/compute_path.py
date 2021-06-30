import numpy as np

def mintree(ngeoms: int,
            geoms: np.ndarray):
  """Computes nearest neighbor tree from origin
  """
  # compute distance between each geometry
  dists = np.zeros((ngeoms,ngeoms))
  for i in range(ngeoms):
    R = geoms[i]
    for j in range(i,ngeoms):
      if i!=j:
        Rnew = geoms[j]
        dists[i,j] = np.sum((R-Rnew)**2.)
        dists[j,i] = dists[i,j]
      else:
        dists[i,i] = np.inf

  # compute nearest neighbors closer to origin
  nns = np.zeros(ngeoms,dtype=int)
  nns[0] = 0
  for i in range(1,ngeoms):
    Ri = dists[i,0]
    dist = dists[i,:]
    idx = np.argsort(dist)
    for j in range(len(idx)):
      Rj = dists[idx[j],0]
      if Rj < Ri:
        nns[i] = idx[j]
        break

  return nns

#def mult_mintree(ngeoms: int,
#                 norigins: int,
#                 geoms: np.ndarray):
#  """Computes nearest neighbor tree from multiple origins
#  """
#  # compute distance between each geometry
#  dists = np.zeros((ngeoms,ngeoms))
#  for i in range(ngeoms):
#    R = geoms[i]
#    for j in range(i,ngeoms):
#      if i!=j:
#        Rnew = geoms[j]
#        dists[i,j] = np.sum((R-Rnew)**2.)
#        dists[j,i] = dists[i,j]
#      else:
#        dists[i,i] = np.inf
#
#  # initialize array with shortest distance paths
#  nns = np.zeros(ngeoms,dtype=int)
#  # "zero" out the origins
#  for i in range(norigins):
#    nns[i] = -9999
#
#  # compute nearest neighbors closer to origin
#  for i in range(norigins,ngeoms):
#    # find closest origin
#    idx0 = np.argsort(dists[i,:norigins])
#    R = dists[i,i]
#    Ri = dists[i,idx0[0]]
#    # get distance of all geometries to this geometry
#    dist = dists[i,:]
#    # get distance of all geometries to this closest origin
#    dist0 = dists[idx0[0],:]
#    # sort them
#    idx00 = np.argsort(dist0)
#    if idx00[0] == i:
#      # current geometry is the closest to this origin
#      nns[i] = idx0[0]
#    else:
#      for j in range(len(idx00)):
#        Rj = dist0[idx00[j]]
#        if Rj < Ri:
#          if dist[idx00[j]] < R:
#            nns[i] = idx00[j]
#            R = dist[idx00[j]]
#
#  return nns
