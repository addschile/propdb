import numpy as np
from propdb import diabatize

if __name__ == "__main__":

  nel  = 2
  nmodes = 2
  modetypes = ['gaussian','gaussian']
  ngeoms = 813
  
  # get adiabatic energies
  eads = np.zeros((ngeoms,2))
  f = open('adpes.dba','r')
  for i in range(ngeoms):
    f.readline() # record no
    es = f.readline().split()
    eads[i,0] = float(es[0])
    eads[i,1] = float(es[1])
  f.close()
  
  # get mode positions
  modes = np.zeros((nmodes,ngeoms))
  f = open('nmode.dba','r')
  for i in range(ngeoms):
    f.readline() # record no
    f.readline() # skip this line
    modes[0,i]  = float(f.readline().split()[0])
    f.readline() # skip this line
    modes[1,i] = float(f.readline().split()[1])
    f.readline() # skip this line
  f.close()
  
  # get nonadiabatic couplings
  f = open('nact.dba','r')
  Fs = np.zeros(ngeoms, dtype=np.ndarray)
  for i in range(ngeoms):
    F = np.zeros((8,3))
    f.readline()
    for j in range(8):
      line = f.readline().split()
      F[j,0] = float(line[0])
      F[j,1] = float(line[1])
      F[j,2] = float(line[2])
    Fs[i] = F
  f.close()
  
  # get geometries
  f = open('geo.dba','r')
  geoms = np.zeros(ngeoms, dtype=np.ndarray)
  for i in range(ngeoms):
    geom = np.zeros((8,3))
    f.readline()
    for j in range(8):
      line = f.readline().split()
      geom[j,0] = float(line[0])
      geom[j,1] = float(line[1])
      geom[j,2] = float(line[2])
    geoms[i] = geom
  f.close()
  
  diabatize(nel,ngeoms,geoms,eads,Fs,krrflag=True,krrargs=[nmodes,modes])

#  diabatize(nel,nmodes,ngeoms,modes,
#  # do krr fit of adiabatic energies
#  alpha = 0.5
#  gam = 1.e-4
#  Kmat = make_gweights_mat(nmodes,ngeoms,alpha,modes,modetypes,gam=gam)
#  Kinv = np.linalg.inv(Kmat)
#  krr_weights = np.zeros((ngeoms,nel))
#  krr_weights[:,0] = np.dot(Kinv,eads[:,0])
#  krr_weights[:,1] = np.dot(Kinv,eads[:,1])
#  vkrr = make_vkrr(nmodes,ngeoms,alpha,modes,modetypes,krr_weights)
#
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
#  # compute nearest neighbors closer to origin
#  nns = np.zeros(ngeoms,dtype=int)
#  nns[0] = -9999
#  for i in range(1,ngeoms):
#    dist = dists[i,:]
#    idx = np.argsort(dist)
#    for j in range(len(idx)):
#      Ri = dists[i,0]
#      Rj = dists[idx[j],0]
#      if Rj < Ri:
#        nns[i] = idx[j]
#        break
#  
#  # do diabatization
#  Cmats = np.zeros((ngeoms,2,2))
#  Cmats[0] = np.eye(2)
#  nsteps = 20
#  ndone = 1
#  for i in range(1,ngeoms):
#    compute_cmat(nns[i],i,nel,nmodes,nns,Cmats,Fs,eads,geoms,modes,modetypes,alpha,krr_weights)
#  
#  # write information to files
#  f = open('krr_opt_db.out','w')
#  fmat = open('krr_opt_trans.out','w')
#  for i in range(ngeoms):
#    # transform adiabatic potential to diabatic
#    vmat = np.zeros((2,2))
#    vmat[0,0] = eads[i,0]
#    vmat[1,1] = eads[i,1]
#    vmat = np.dot(Cmats[i].T, np.dot(vmat, Cmats[i]))
#    # write diabatic energies
#    f.write('#Record No. : %d\n'%(i+1))
#    f.write('%.8f '%(vmat[0,0]))
#    f.write('%.8f '%(vmat[0,1]))
#    f.write('%.8f\n'%(vmat[1,1]))
#    f.flush()
#    # write transformation matrix
#    fmat.write('#Record No. : %d\n'%(i+1))
#    fmat.write('%.8f %.8f\n'%(Cmats[i,0,0],Cmats[i,0,1]))
#    fmat.write('%.8f %.8f\n'%(Cmats[i,1,0],Cmats[i,1,1]))
#    fmat.flush()
#  f.close()
#  fmat.close()
#
#  f = open('krr_opt_db.out','r')
#  for i in range(ngeo):
#    f.readline() # record no
#    es = f.readline().split()
#    ead1_krr[i]  = float(es[0])
#    ead12_krr[i] = float(es[1])
#    ead2_krr[i]  = float(es[2])
#  f.close()
#  
#  # plot diabatized pes
#  fig = plt.figure()
#  ax = fig.add_subplot(111, projection='3d')
#  ax.scatter(q5[:ngeo],q14[:ngeo],ead1_krr,'b')
#  ax.scatter(q5[:ngeo],q14[:ngeo],ead2_krr,'r')
#  ax.set_xlabel('q5')
#  ax.set_ylabel('q14')
#  plt.savefig('pes.png')
#  plt.show()
#
#  # plot diabatic coupling
#  fig = plt.figure()
#  ax = fig.add_subplot(111, projection='3d')
#  ax.scatter(q5[:ngeo],q14[:ngeo],ead12_krr,'b')
#  ax.set_xlabel('q5')
#  ax.set_ylabel('q14')
#  plt.savefig('db_coup.png')
#  plt.show()
