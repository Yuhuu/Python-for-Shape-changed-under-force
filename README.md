
<h2>Here are the tasks </h2>
Using pen & paper, add element stiffness matrix KE to global system
stiffness matrix KS at global degrees of freedom (DOFs) given by L2G.

KE = array([[1 ,2 ,3 ,4 ],
            [5 ,6 ,7 ,8 ],
            [9 ,10,11,12],
            [13,14,15,16]])*1.

KS = array([[0,0,0,0,0,1],
            [1,1,1,1,0,0],
            [1,1,1,1,0,1],
            [1,1,1,1,0,0],
            [1,1,1,1,1,0],
            [1,1,1,1,0,1]])*1.

L2G = array( [ 5, 6, 3, 4 ] )


Function header:
def add2stiffness(ke,local2global,ks):
    # ke           : 4 x 4  element matrix : input
    # local2global : vector with length 4, pointing from local to global DOF : input
    # example local2global = array( [ 5, 6, 3, 4 ] )
    # ks           : system stiffness matrix : input and output
    # loop over columns
      # loop over rows
        # add element stiffness component to system stiffness component
    return ks
    
    Use pen & paper and modify the equation system (1) for one displacement
boundary condition with one prescribed displacement "Didof" for degree of
freedom (DOF) number "idof".

  (1)         [KS]{D} = {F}

      [KS] is the system stiffness matrix.
      {F}  is the right-hand-side, where external forces already have been applied.

KS = array([[1.,2.,3.],
            [4.,5.,6.],
            [7.,8.,9.]])   # System matrix
F=array([11,12,13])*1.     # Applied forces, column vector
idof = 2                   # Known DOF
Didof = 3.                 # Known displacement at known DOF
#-------------------------------------------------------------------------------


Write a Python function modifying an equation system like (1) for one
displacement boundary condition with one prescribed displacement Didof for 
degree of freedom number idof, and show use of the function with the inputs
KS,F,idof,Didof given in problem 2.1) to verify that it works properly!

#Function header:
def enforce_displacement(ks,F,idof,Didof):
  # check dimensions
  # move terms involving Didof to the RHS:
   # - subtract terms involving Didof from RHS
   # - blank column idof-1 in ks
   # - blank row idof-1 in ks
  # set diagonal element of ks to 1.0
  # set RHS row idof-1 identical to prescribed displacement Didof

  # WARNING: Python indexing of arrays and matrices start with 0 !

  return [ks,F]
