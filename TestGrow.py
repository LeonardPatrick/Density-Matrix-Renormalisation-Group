import numpy as np
from scipy.sparse import kron
import time
from scipy.linalg import eigh    # importing dense eigensolver
from scipy.sparse.linalg import eigsh # importing Arpack
from numpy import * # to include zeros function



# h is size of submatrix: hxh.
# rho is the full density matrix.
# take the example of A[0:3, 0:3]. The first argument is for the rows; 0:3 means that you start at row 0 (the first row) and you go down to the third row (index row 2), so that one is dealing with three rows. Similarly, the second argument corresponds with columns and 0:3 goes form column index = 0 to column index = 2. 

#RDM = np.array()


def RDMfunc(h, rho):

    # h =2   # the reduced density matrix will be  2x2. h also represents the dimension of the Hilbert space of the system. In this case h = 2 indicates that we are dealing with a single site of spin-1/2\. Since the matrix is 8x8, it must mean that the environment is represented by 4x4 matrix, and thus corresponds to two sites.                                                                             

     n = len(rho)/h  # this represents the dimension of the submatrices used to create the reduced density matrix                                                                                         
     global RDM # global because you wish to make the array RDM available outside the function                                                                                                               
     RDM = zeros((h,h))

     for i in range (0,h):
       for j in range (0,h):
        RDM[i,j] = np.trace( rho[i*n:i*n + n, j*n:j*n + n])




         



start_time = time.time()



# Construct Sz, S+ and S- for a single site 

Sz = np.array([[ 0.5, 0.0], [0.0 ,-0.5]], dtype = 'd')

Sp = np.array([[0.0, 1.0], [0.0, 0.0] ], dtype = 'd')

Sm = np.array([ [ 0.0, 0.0],[1.0, 0.0]], dtype = 'd')





#Initial Hamiltonians for left and right sides; two site Hamiltonians.

#HL = zeros((4,4))

HL =  np.kron(Sz, Sz) + 0.5 * (np.kron(Sp, Sm) + np.kron(Sm, Sp) )

HR = np.kron(Sz, Sz) + 0.5 * (np.kron(Sp, Sm) + np.kron(Sm, Sp) )

                                                 


SztL = Sz;
SptL = Sp;
SmtL = Sm;


SztR = Sz;
SptR = Sp;
SmtR = Sm;

  

# Growing left and right blocks to a desired size inside for loop; outside of for loop, superblock is created by connecting left and right blocks to two centre sites. 


# Recursion 


b = 5

N = b + 1

LN = 2.0 * (N - 1.0) 




 



for i in range(3, N):
 
 # left block
      SztL = np.kron( np.identity(2), SztL)
      SptL = np.kron( np.identity(2), SptL)
      SmtL = np.kron( np.identity(2), SmtL)



 # right block
      SztR = np.kron(SztR, np.identity(2))
      SptR = np.kron(SptR, np.identity(2))
      SmtR = np.kron(SmtR, np.identity(2))




 # Refer to p. 1566 of notes

      HL = np.kron(HL, np.identity(2)) + np.kron(SztL, Sz) + 0.5 *( np.kron(SptL, Sm) + np.kron(SmtL, Sp)) 

      HR = np.kron(Sz, SztR) + 0.5 * (np.kron(Sp, SmtR) + np.kron(Sm, SptR)) + np.kron(np.identity(2), HR )








DL = 2**(N - 2)

DR = 2**(N - 2)




# Left

SzSHL = np.kron(identity(DL), Sz)
SpSHL = np.kron(identity(DL), Sp)
SmSHL = np.kron(identity(DL), Sm)





# Right



SzSHR = np.kron(Sz, identity(DR))
SpSHR = np.kron(Sp, identity(DR))
SmSHR = np.kron(Sm, identity(DR))



H = np.kron(HL, identity(DR * 2)) + np.kron(identity( DL * 2), HR ) + np.kron(SzSHL, SzSHR) + 0.5 * (np.kron(SpSHL, SmSHR)) + 0.5 * (np.kron(SmSHL, SpSHR))



print np.shape(H)

print ' '

GSE, GS = eigsh(H, k = 1, which = "SA"  )   



print GSE/(b*2)







# Forming superblock. Refer to p. 1568 of notes

DL = 2**(N - 2)

DR = 2**(N - 2)


# From p. 1570 in in notes



# Left

SzSHL = np.kron(identity(DL), Sz)
SpSHL = np.kron(identity(DL), Sp)
SmSHL = np.kron(identity(DL), Sm)





# Right



SzSHR = np.kron(Sz, identity(DL))
SpSHR = np.kron(Sp, identity(DL))
SmSHR = np.kron(Sm, identity(DL))



# Final Hamiltonian/ Super Hamiltonian (refer to p. 1568 of notes)

# Refer to Feiguin's notes p. 49  (19) as to why the dimensions of the intial H are (DR * DL * 4, DR * DL * 2).



H = np.kron(HL, identity(DR * 2)) + np.kron(identity( DL * 2), HR ) + np.kron(SzSHL, SzSHR) + 0.5 * (np.kron(SpSHL, SmSHR)) + 0.5 * (np.kron(SmSHL, SpSHR))



print np.shape(H)

print " "

# Find Ground state



GSE, GS = eigh(H )






print " "

print LN
print " "
print GSE

print " "


