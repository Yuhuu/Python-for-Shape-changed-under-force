
 from pylab import *
import numpy as np
np.set_printoptions(precision=5,suppress=True)

def stiffness_3D_rod(X1,Y1,Z1,X2,Y2,Z2,E,A):
  L = ((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)**0.5
  s = (X2-X1)/L
  m = (Y2-Y1)/L
  n = (Z2-Z1)/L
  K=array([[ s**2, s*m , s*n ,-s**2,-s*m ,-s*n ],
           [ s*m , m**2, m*n ,-s*m ,-m**2,-m*n ],
           [ s*n , m*n , n**2,-s*n ,-m*n ,-n**2],
           [-s**2,-s*m ,-s*n , s**2, s*m , s*n ],
           [-s*m ,-m**2,-m*n , s*m , m**2, m*n ],
           [-s*n ,-m*n ,-n**2, s*n , m*n , n**2]])*E*A/L
  return K

def add2system(ke,l2g,ks):
    for i in range(len(l2g)):
        for j in range(len(l2g)):
            ks[(l2g[i]-1),(l2g[j]-1)]=ks[l2g[i]-1,l2g[j]-1]+ke[i,j]
    return ks

def enforce_displacement(ks,F,idof,didof):
    for i in range(len(idof)):
        F[:]=F[:]-ks[:,idof[i]-1]*didof[i]
        ks[:,idof[i]-1]=0
    for j in range(len(idof)):
        ks[idof[j]-1,:]=0
    for l in range(len(idof)):
        ks[idof[l]-1,idof[l]-1]=1
        F[idof[l]-1]=didof[l]
    return ks,F

coordinates = array ([[0,-1, 0],      # x,y Coordinates for node 1
                      [0, 1, 0],
                      [2,-1, 0],
                      [2, 1, 0],
                      [4,-1, 0],
                      [4, 1, 0],
                      [6,-1, 0],
                      [6, 1, 0],
                      [8, 0, 0],
                      [6, 0, 1.5],
                      [5, 0, 1.5],
                      [3, 0, 1.5],
                      [1, 0, 1.5],
                      [0, 0, 1.5]])     # x,y Coordinates for node 3

l2ge = array ([[1,2],                   # Connecting nodes in element 1,2,3(node to
node)
               [1,3],
               [3,4],
               [3,5],
               [5,6],
               [5,7],
               [7,8],
               [7,9],
               [2,4],
               [4,6],
               [6,8],
               [8,11],
               [7,11],
               [6,11],
               [5,11],
               [6,12],
               [5,12],
               [4,12],
               [3,12],
               [4,13],
               [3,13],
               [2,13],
               [1,13],
               [9,10],
               [8,10],
               [7,10],
               [1,13],
               [10,11],
               [11,12],
               [12,13],
               [13,14],
               [2,3],
               [3,6],
               [6,7],
               [8,9]])
#              b[mm] h[mm] E[N/mm^2] rho[kg/m^3]
#Mat = array ([[ 210000.,100.], # b,h, E-modulus and rho for element 1,2,3
#              [  70000.,120.],
#              [  50000.,200.]])

E = 210000
A = 800

nelements = size(l2ge,axis=0)           # Calculating number of elements
nnodes = size(coordinates,axis=0)       # Calculating number of nodes

ndof = 3                                # Degrees of freedom per node
nn_ele = 2                              # Numer of nodes per element

idof = [1,2,3,4,5,6,40,41,42]                        # Known degrees of freedom for
displacements
Didof = [0.,0.,0.,0.,0.,0.,0.,0.,0]                  # Magnitude of the known
displacements

f_idof = [27]                          # Known degrees of freedom for forces
Fidof = [-5000*9.81]                  # Magnitude of the known forces
F = zeros(ndof*nnodes)                  # Creating the force matrix

for x in range (size(Fidof)):           # This loop can also be created as a function
    F[f_idof[x]-1] = Fidof[x]
    KEG = zeros([ndof*nn_ele,ndof*nn_ele,nelements])    # Zero matrises for the
kegs, stacking in 3D
    KS = zeros([ndof*nnodes,ndof*nnodes])   # Global matrix
    l2g = zeros([ndof*nn_ele,nelements])    # Local 2 global matrix
for x in range (nelements):             # Creating all KEGs and KS with one loop
    x1 = coordinates[l2ge[x,0]-1,0]     # Picking data for the element x according
to the l2g matrix
    y1 = coordinates[l2ge[x,0]-1,1]
    z1 = coordinates[l2ge[x,0]-1,2]
    x2 = coordinates[l2ge[x,1]-1,0]
    y2 = coordinates[l2ge[x,1]-1,1]
    z2 = coordinates[l2ge[x,1]-1,2]
#    E = Mat[x,0]
#    A = Mat[x,1]
    KEG[:,:,x] = stiffness_3D_rod(x1,y1,z1,x2,y2,z2,E,A) # Creating the local
stiffness matrix for element x
    l2g[0,x]= l2ge[x,0]*ndof-2                      # Expanding Local 2 global from
nodes to dofs
    l2g[1,x]= l2ge[x,0]*ndof-1
    l2g[2,x] =l2ge[x,0]*ndof
    l2g[3,x]= l2ge[x,1]*ndof-2
    l2g[4,x] =l2ge[x,1]*ndof-1
    l2g[5,x]= l2ge[x,1]*ndof
    KS = add2system(KEG[:,:,x],l2g[:,x],KS)     # adding the KEG to KS
KSorg = copy(KS)                                    # Stores the original KS before
tampering with enforce displacement
[KS,F] = enforce_displacement(KS,F,idof,Didof)  # Manipulates KS with respect to
idof and Didof to make the system solvable
KSinv = inv(KS)
d = dot(KSinv,F)                                    # Computing displacements
F_reaction = dot(KSorg,d)                           # Computing the reaction forces

print '\nDisplacements [mm]:'
print d
print '\nReaction forces [N]:'
print F_reaction

RF_1_3 = dot(KEG[:,:,2],d[[0,1,2,6,7,8]])

RF_3_5 = dot(KEG[:,:,4],d[[6,7,8,12,13,14]])

RF_12_13 = dot(KEG[:,:,29],d[[33,34,35,36,37,38]])

print 'Axial element forces between keypoiunt 1 and 3', RF_1_3, '[N] \n'