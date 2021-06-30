import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

q5  = np.zeros(813)
q14 = np.zeros(813)
f = open('nmode.dba','r')
for i in range(813):
    f.readline() # record no
    f.readline() # skip this line
    q5[i]  = float(f.readline().split()[0])
    f.readline() # skip this line
    q14[i] = float(f.readline().split()[1])
    f.readline() # skip this line
f.close()

ead1  = np.zeros(813)
ead2  = np.zeros(813)
ead12 = np.zeros(813)
f = open('pes.dba','r')
for i in range(813):
    f.readline() # record no
    es = f.readline().split()
    ead1[i]  = float(es[0])
    ead12[i] = float(es[1])
    ead2[i]  = float(es[2])
f.close()

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(q5,q14,ead1,'b')
#ax.scatter(q5,q14,ead2,'r')
#ax.set_xlabel('q5')
#ax.set_ylabel('q14')
##plt.savefig('pes_og.png')
#plt.show()
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(q5,q14,ead12,'b')
#ax.set_xlabel('q5')
#ax.set_ylabel('q14')
##plt.savefig('db_coup_og.png')
#plt.show()

ead1  = np.zeros(813)
ead2  = np.zeros(813)
ead12 = np.zeros(813)
f = open('db.out','r')
for i in range(813):
    f.readline() # record no
    es = f.readline().split()
    ead1[i]  = float(es[0])
    ead12[i] = float(es[1])
    ead2[i]  = float(es[2])
f.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(q5,q14,ead1,'b')
ax.scatter(q5,q14,ead2,'r')
ax.set_xlabel('q5')
ax.set_ylabel('q14')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(q5,q14,ead12,'b')
ax.set_xlabel('q5')
ax.set_ylabel('q14')
plt.show()
