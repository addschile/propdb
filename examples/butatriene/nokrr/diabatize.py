import numpy as np
from propdb import diabatize

if __name__ == "__main__":

  nel  = 2
  nmodes = 2
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
    Fs[i] = F/(eads[i,1]-eads[i,0])
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
  
  diabatize(nel,ngeoms,geoms,eads,Fs)
