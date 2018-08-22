from __future__ import print_function, division  # requires Python >= 2.6 

import numpy as np
from scipy.sparse import kron, identity
from scipy.sparse.linalg import eigsh  # Lanczos routine from ARPACK                                                                                                                                     
import time

np.set_printoptions(precision=10, suppress=True, threshold=10000, linewidth=300)

start_time = time.time()

         

def truncEr(array, n):

     sum = 0.0

     for i in range(0,n):
          sum = sum + array[i]
     
     return 1 - sum    




'''
print(" ")
print(" Trunc")
print( truncEr(evals, m))
print(" ")
'''


# list/array to store matrices for finite algorithm; refer to p. 1691 of notes.

HLMAT = []
HRMAT = [] 



# Construct Sz, S+ and S- for a single site 

Sz = np.array([[ 0.5, 0], [0 ,-0.5]], dtype = 'd')
Sp = np.array([[0, 1], [0, 0] ], dtype = 'd')
Sm = np.array([ [ 0, 0],[1, 0]], dtype = 'd')



#Initial Hamiltonians for left and right sides; two site Hamiltonians.

HL = kron(Sz, Sz) + 0.5 * (kron(Sp, Sm) + kron(Sm, Sp) )
HR = kron(Sz, Sz) + 0.5 * (kron(Sp, Sm) + kron(Sm, Sp) )


# Adding HR to the list; refer to p. 1691 of notes  

HRMAT.append(HR)
HLMAT.append(HL)


# Initialising matrices for growth

SztL = Sz;
SptL = Sp;
SmtL = Sm;


SztR = Sz;
SptR = Sp;
SmtR = Sm;


  




###############################################
# Forming superblock. Refer to p. 1568 of notes


m = 4 # How many states to keep associated with reduced density matrix of two site system. In this case, all states are kept.
LN = 4 # Total number of lattice sites.


DL = 2
DR = 2


# From p. 1570 in in notes

# Left

SzSHL = kron(identity(DL), Sz)
SpSHL = kron(identity(DL), Sp)
SmSHL = kron(identity(DL), Sm)


# Right

SzSHR = kron(identity(DR), Sz)
SpSHR = kron(identity(DR), Sp)
SmSHR = kron(identity(DR), Sm)




# Final Hamiltonian/ Super Hamiltonian (refer to p. 1568 of notes)

H = kron(HL, identity(DR * 2)) + kron(identity(DL * 2), HR) +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR) #works



# Find Ground state

(GSE,), GS = eigsh(H, k=1, which="SA")


print (LN)
print (GSE/LN)
print (" ")




# Construct reduced density matrix

psi = GS.reshape([2**(LN/2), -1], order='C')
rho = np.dot(psi, psi.conjugate().transpose())



# Diagonalise reduced density matrix

evals, evecs = np.linalg.eigh(rho)




# Form the truncated reduced density matrix


idx = evals.argsort()  # Order eigenvalues 
idx = idx[::-1]        # Reverse order so that largest eigenvalues are first
U = evecs[:,idx[:m]]   # Create transformation matrix with m^th largest eigenvalues  


 

# Initialsing spin operators that connect two blocks for growth
     
SztL = kron(identity( 2 ) , SztL)
SptL = kron(identity( 2 ) , SptL)
SmtL = kron(identity( 2 ) , SmtL)

          
# right block. Refer to. 1676
     

SztR = kron(SztR,identity(2))
SptR = kron(SptR,identity(2))
SmtR = kron(SmtR,identity(2))




M  = 10

N = 3

L = N + 8 # The number 8 indicates how many elements will be in the list/ array HRMAT due to the following loop. However, in total, there will be 9 since we have the contribution of HR above with the lattice possessing 4 sites. 

m = 4  # m #initialising ml to be that of m at the start of the code.




for p in range(N,L):      
      


      DL = m # new dimension of truncated block. ml will vary depending on whether the dimension of the Hilbert space is smaller or larger than a fixed size M (see below)
      DR = m
      LN = 2.0 * (p)     #and not 2*(p - 1) as in the part of the algorithm that creates the intial lattice.  


             
      # transformation of left block
  
      SztL = U.conjugate().transpose().dot(SztL.dot(U)) 
      SptL = U.conjugate().transpose().dot(SptL.dot(U))
      SmtL = U.conjugate().transpose().dot(SmtL.dot(U))


       

      # transformation of right block

      SztR = U.conjugate().transpose().dot(SztR.dot(U))
      SptR = U.conjugate().transpose().dot(SptR.dot(U))
      SmtR = U.conjugate().transpose().dot(SmtR.dot(U))


      HL = U.conjugate().transpose().dot(HL.dot(U))
      HR = U.conjugate().transpose().dot(HR.dot(U))
 

     

      HL = kron(HL, identity(2)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 
      HR = kron(HR, identity(2)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 



      HRMAT.append(HR) # Adding HR to the list; refer to p. 1691 and 1708 of notes. 
      HLMAT.append(HL) # Adding HR to the list; refer to p. 1691 and 1708 of notes  



      #create superblock

      # left
  

      SzSHL = kron(identity(DL), Sz)  # Appear to be correct
      SpSHL = kron(identity(DL), Sp)
      SmSHL = kron(identity(DL), Sm)

     


      # right


      
      SzSHR = kron(identity(DR), Sz)  # Appear to be correct
      SpSHR = kron(identity(DR), Sp)
      SmSHR = kron(identity(DR), Sm)
      



      H = kron(HL, identity(DR * 2))  + kron(identity(DL * 2), HR) +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR) #works





      # Find Ground state


      (GSE,), GS = eigsh(H, k=1, which="SA")


      print (" ")
      print ("LN = ", LN )
      print ("E/L = ", GSE/LN)
      print (" ")




      # Construct reduced density matrix

      if 2**(p-1) < M:
           sd = 2**(LN/2)

      else:
           sd = 2*M    


      psi = GS.reshape([sd, -1], order='C')    # Problem here when truncaiton begins. If statement needs to be implemented in order to ensure that the wavefunction is reshaped correctly.
      RDM = np.dot(psi, psi.conjugate().transpose())
      evalsl, evecsl = np.linalg.eigh(RDM)


      # Form the truncated reduced density matrix

      idx = evalsl.argsort()
      idx = idx[::-1] 


      if 2**p < M:
           m = 2**p
      else:
           m = M

      U = evecsl[:,idx[:m]]  # This works

          




      # Growth of operators for next iteration

      if 2**p < M:
 
           SztL = kron(identity( 2**(p-1) ) , Sz)
           SptL = kron(identity( 2**(p-1) ) , Sp)
           SmtL = kron(identity( 2**(p-1) ) , Sm)

           SztR = kron(identity( 2**(p-1) ) , Sz)
           SptR = kron(identity( 2**(p-1) ) , Sp)
           SmtR = kron(identity( 2**(p-1) ) , Sm)

      else:

           # Grow SztL, SmtR, etc. Refer to p. 1676

           # left
  

           SztL = kron(identity(DL), Sz) 
           SptL = kron(identity(DL), Sp)
           SmtL = kron(identity(DL), Sm)
          
           # right 
         
           SztR = kron(identity(DR), Sz)
           SptR = kron(identity(DR), Sp) 
           SmtR = kron(identity(DR), Sm)


      




##############################################################################

# Finite part

# Use HL from the last part of algorithm to grow left side.
# If you are ending with L = 20, it means that the system size is 10 sites. The following then will grow the left side by one site to make the system 11 sites. Therefore the environment is 9 sites. In turn this means that one should return to the case of L = 18 and store the right hand side block.   


# Start with symmetric configuration from infinite algorithm ---------00---------
# Grow left system by one site and shrink right system by one site. Keep in mind that the length of the chain remains at 20 sites.



# The following will end with ------------------00

for i in range(1,9):    # loop starts at 1 rather than 0 because the indexing of HRMAT is -i-1 and the first element to be accessed should be HR9 (since we have 20 sites in total [10 initially on both sides], left grows by one and right shrinks by one, so it has 9 sites [refer to p. 1691 of notes]).

    

    print(" ")  
    print("Number of sites on the right = ", 10 - i )
    print(" ")



    if i == 6:
        DR = 8

    if i == 7:
        DR = 4

    if i == 8:
        DR = 2  # the last change of DRm before the sweep returns in the opposite direction.



    # Transformation of spin operators connecting left block to right block. 

    SztL = U.conjugate().transpose().dot(SztL.dot(U))  
    SptL = U.conjugate().transpose().dot(SptL.dot(U))
    SmtL = U.conjugate().transpose().dot(SmtL.dot(U))




    # Transformation of left system block

    HL = U.conjugate().transpose().dot(HL.dot(U))



    # Growth of left side by one site
           
    HL = kron(HL, identity(2)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 




   # HRMAT.append(HR) # Adding HR to the list; refer to p. 1691 and 1708 of notes. 
    HLMAT.append(HL)





    if i == 7:
        HLNEW = HL 



    # left
  
    SzSHL = kron(identity(DL), Sz)
    SpSHL = kron(identity(DL), Sp)
    SmSHL = kron(identity(DL), Sm)

     

    # right
      
    SzSHR = kron(identity(DR), Sz)
    SpSHR = kron(identity(DR), Sp)
    SmSHR = kron(identity(DR), Sm)




    # i == 8 signifies that the right system is now only composed of two sites. At this stage the right side becomes the system. In this current form of the code, this means that the superblock Hamilto      nian must be rearranged. Refer to p. 1700

    if i == 8:

        H = kron(HRMAT[-i-1], identity(DL * 2)) + kron(identity(DR * 2), HL)  +  0.5 * (kron(SpSHR, SmSHL)) +  0.5 * (kron(SmSHR, SpSHL)) + kron(SzSHR, SzSHL)

    else:
  
        H = kron(HL, identity(DR * 2)) + kron(identity(DL * 2), HRMAT[-i-1])  +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR) #works

       
   

    (GSE,), GS = eigsh(H, k=1, which="SA")




    print (" ")
    print ("LN = ", LN )
    print ("E/L = ", GSE/20) # 20 is used because the lattice only possesses 20 sites.
    print (" ")

      

    if i ==8:

        sd = 4
   
   

    # Construct reduced density matrix

    psi = GS.reshape([sd, -1], order='C')      
    RDM = np.dot(psi, psi.conjugate().transpose())
    evalsl, evecsl = np.linalg.eigh(RDM)

    idx = evalsl.argsort()
    idx = idx[::-1] 



    m = M

    if i == 8:
        m = 4



    U = evecsl[:,idx[:m]]  


    # left
  

    SztL = kron(identity(DL), Sz)
    SptL = kron(identity(DL), Sp)
    SmtL = kron(identity(DL), Sm)   
    

    # right    
      
    SztR = kron(identity(DR), Sz)
    SptR = kron(identity(DR), Sp) 
    SmtR = kron(identity(DR), Sm)



















# The following will end with 00------------------


#######################################
# Sweeping in the opposite direction


HR = HRMAT[0] # Initial conidition
HRMAT = []  # The list is emptied as previous matrix are no longer required, except for the first element representing a two-site subsystem; this is resolved in the next line.  
HRMAT.append(HR)


M = 10

for i in range(1,17):   



    print(" ")  
    print("Number of sites on the right = ", i + 2 )
    print(" ")




    if i == 1:
        DR = 4



    if i == 2:
        DR = 8
    


    if i == 3:
        DR = 10



    if i > 3:
        DR = 10


    if i == 14:
        DL = 8


    if i == 15:
        DL = 4



    if i == 16:
        DL = 2






    # Transformation of right side operators  

    SztR = U.conjugate().transpose().dot(SztR.dot(U)) 
    SptR = U.conjugate().transpose().dot(SptR.dot(U))
    SmtR = U.conjugate().transpose().dot(SmtR.dot(U))





    HR = U.conjugate().transpose().dot(HR.dot(U))   # HR instead of HL

    # Growth
           
    HR = kron(HR, identity(2)) + kron(SztR, Sz) + 0.5 *( kron(SptR, Sm) + kron(SmtR, Sp)) 


    HRMAT.append(HR)






    # left
  
    SzSHL = kron(identity(DL), Sz) 
    SpSHL = kron(identity(DL), Sp)
    SmSHL = kron(identity(DL), Sm)

     

    # right
      
    SzSHR = kron(identity(DR), Sz) 
    SpSHR = kron(identity(DR), Sp)
    SmSHR = kron(identity(DR), Sm)






    if i == 16:

        H = kron(HLMAT[-i-1], identity(DR * 2)) + kron(identity(DL * 2), HR)  +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR) # Strange: this seems to work. Note now that the HRMAT is swapped with HL in the special case that there are two sites remaining on the right (i == 8)  . Note that all other matrices are swapped, and DLm and DRm (compare with below). What does this mean? REFER TO P. 1700.  


    else:
  
        H = kron(HR, identity(DL * 2)) + kron(identity(DR * 2), HLMAT[-i-1])  +  0.5 * (kron(SpSHR, SmSHL)) +  0.5 * (kron(SmSHR, SpSHL)) + kron(SzSHR, SzSHL) # Strange: this seems to work. Note now that the HRMAT is swapped with HL in the special case that there are two sites remaining on the right (i == 8)  . Note that all other matrices are swapped, and DLm and DRm (compare with below). What does this mean? REFER TO P. 1700.  



    # Find ground state       

    (GSE,), GS = eigsh(H, k=1, which="SA")


    print(" ")
    print(GSE/20)
    print(" ")
    


    # Construct reduced density matrix

    if 2*(DR - 4) < M:
       sd = 2 * DR
    elif i ==16:
        sd = 4
    else:
       sd = 2*M




 
    psi = GS.reshape([sd, -1], order='C')    # Problem here when truncaiton begins. If statement needs to be implemented in order to ensure that the wavefunction is reshaped correctly.
    RDM = np.dot(psi, psi.conjugate().transpose())

    evalsl, evecsl = np.linalg.eigh(RDM)

    idx = evalsl.argsort()
    idx = idx[::-1] 


    # Form truncated matrix

    if 2*(DR) < M:  # No -4, as above with rs, since 16 > MML and so we wish to truncate.
       m = 2 * DR 
    else:
       m = M

    

    U = evecsl[:,idx[:m]] 


    # left
  
    SztL = kron(identity(DL), Sz) 
    SptL = kron(identity(DL), Sp)
    SmtL = kron(identity(DL), Sm)


    # right block. Refer to. 1676
             
    SztR = kron(identity(DR), Sz)
    SptR = kron(identity(DR), Sp) 
    SmtR = kron(identity(DR), Sm)


    


 
 
 



###################################################################### 
# Sweeping back in the original direction
# Now the left side becomes the SYSTEM



# The following will end with ---------00---------


# Transformation matrix in the last part of previous section must be transformed as to be 4x4 to reflect that the system is now represented by the left side (which is 2 sites for the last part).







HL = HLMAT[0]  # Initialising
HLMAT = []
HLMAT.append(HL)

DR = 10
DL = 4



for i in range(1,9):  

   





    if i == 1:
        DL = 4

    if i == 2:
        DL = 8

    if i == 3:
        DL =10 


   
   

    SztL = U.conjugate().transpose().dot(SztL.dot(U)) 
    SptL = U.conjugate().transpose().dot(SptL.dot(U))
    SmtL = U.conjugate().transpose().dot(SmtL.dot(U))



    HL = U.conjugate().transpose().dot(HL.dot(U))   # HR instead of HL




    # Growth
           
    HL = kron(HL, identity(2)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 

       
   
  

    HLMAT.append(HL)

    # left
  
    SzSHL = kron(identity(DL), Sz) 
    SpSHL = kron(identity(DL), Sp)
    SmSHL = kron(identity(DL), Sm)

     
    # right
    
    SzSHR = kron(identity(DR), Sz) 
    SpSHR = kron(identity(DR), Sp)
    SmSHR = kron(identity(DR), Sm)




  
  

    
    

    #Find HRMAT. HRNEW below is the HR matrix in the previous stage of this algorithm; it is the matrix associated with the penultimate run. 

    H = kron(HL, identity(DR * 2)) + kron(identity(DL * 2), HRMAT[-i-1]) +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR)





    # Find ground state

    (GSE,), GS = eigsh(H, k=1, which="SA")

    print(" ")
    print(GSE/20)
    print(" ")
    





    # Construct reduced density matrix

    if 2*(DL - 4) < M:
       sd = 2 * DL
    else:
       sd = 2*M



 
    psi = GS.reshape([sd, -1], order='C')    # Problem here when truncaiton begins. If statement needs to be implemented in order to ensure that the wavefunction is reshaped correctly.
    RDM = np.dot(psi, psi.conjugate().transpose())

    evalsl, evecsl = np.linalg.eigh(RDM)




   # print(" ")
   # print("shape")
   # print(np.shape(H))
   # print(" ")


    idx = evalsl.argsort()
    idx = idx[::-1] 


    # Form truncated matrix

    if 2*(DL) < M:  # No -4, as above with rs, since 16 > MML and so we wish to truncate.
       m = 2 * DL 
    else:
       m = M


    

    U = evecsl[:,idx[:m]] 



    # left
  
    SztL = kron(identity(DL), Sz) 
    SptL = kron(identity(DL), Sp)
    SmtL = kron(identity(DL), Sm)


    # right block. Refer to. 1676
             
    SztR = kron(identity(DR), Sz)
    SptR = kron(identity(DR), Sp) 
    SmtR = kron(identity(DR), Sm)




print (time.time() - start_time, " seconds")
