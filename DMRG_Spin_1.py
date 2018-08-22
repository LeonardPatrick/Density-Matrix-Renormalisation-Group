from __future__ import print_function, division  # requires Python >= 2.6 

import numpy as np
from scipy.sparse import kron, identity
from scipy.sparse.linalg import eigsh  # Lanczos routine from ARPACK                                                                                                                                     
import math
import time

#np.set_printoptions(precision=12, suppress=True, threshold=10000, linewidth=300)
np.set_printoptions(precision=10, suppress=True, threshold=10000, linewidth=300)
start_time = time.time()

         

def truncEr(array, n):

     sum = 0.0

     for i in range(0,n):
          sum = sum + array[i]
     
     return 1 - sum    








################################################
# INFINITE SYSTEM ALGORITHM


# list/array to store matrices for finite algorithm; refer to p. 1691 of notes.

HLMAT = []
HRMAT = [] 



# list/ array to store transformation/ truncation matrices for wave function transformation.

UENVL = []
UENVR = []




# Construct Sz, S+ and S- for a single site 

sqtwo=1.41421356237309504880168872421

Sz = np.array([[ 1.0, 0, 0], [0 , 0, 0], [0, 0 , -1.0]], dtype = 'd')
Sp = np.array([[0, sqtwo , 0], [0, 0, sqtwo ], [0, 0, 0] ], dtype = 'd')
Sm = np.array([[0, 0, 0], [sqtwo , 0, 0], [0, sqtwo, 0] ], dtype = 'd')








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


dim = 3 # Dimension of Hilbert space associated wiht a single site 



m = dim**2 # How many states to keep associated with reduced density matrix of two site system. In this case, all states are kept.
LN = 4 # Total number of lattice sites.


DL = dim
DR = dim


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

H = kron(HL, identity(DR * dim)) + kron(identity(DL * dim), HR) +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR) #works



# Find Ground state

(GSE,), GS = eigsh(H, k=1, which="SA")






print ("LN = ", LN)
print ("E/L = ",GSE/LN)
print (" ")




# Construct reduced density matrix



sd = HL.shape[0]

psi = GS.reshape([sd, -1], order='C')




rho = np.dot(psi, psi.conjugate().transpose())



# Diagonalise reduced density matrix

evals, evecs = np.linalg.eigh(rho)




# Form the truncated reduced density matrix


idx = evals.argsort()  # Order eigenvalues 
idx = idx[::-1]        # Reverse order so that largest eigenvalues are first
U = evecs[:,idx[:m]]   # Create transformation matrix with m^th largest eigenvalues  



UENVL.append(U)         # Storing matrices for wave function transformation. Refer to p. 1758

 

# Initialsing spin operators that connect two blocks for growth
     
SztL = kron(identity( dim ) , SztL)
SptL = kron(identity( dim ) , SptL)
SmtL = kron(identity( dim ) , SmtL)

          
# right block. Refer to. 1676
     

SztR = kron(SztR,identity(dim))
SptR = kron(SptR,identity(dim))
SmtR = kron(SmtR,identity(dim))



M = 10

N = 3

Sit = 20               # Number of sites required to output at an end of infinite algorithm;  must be an even number

It = (Sit/2) - 2  # No of iterations required to produce number of sits Sit. Refer to p. 1712 
It2 = Sit - (Sit/2) - 2
It3 = Sit - 4 # refer to p. 1715
It4 = (Sit/2) - 2
# Must convert It to integer for it to be used in for loop.

L = N + int(It)

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
 

     

      HL = kron(HL, identity(dim)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 
      HR = kron(HR, identity(dim)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 



      HRMAT.append(HR) # Adding HR to the list; refer to p. 1691 and 1708 of notes. 
      HLMAT.append(HL) # Adding HR to the list; refer to p. 1691 and 1708 of notes  



      #create superblock

      # left
  

      SzSHL = kron(identity(DL), Sz)
      SpSHL = kron(identity(DL), Sp)
      SmSHL = kron(identity(DL), Sm)

     


      # right


      
      SzSHR = kron(identity(DR), Sz) 
      SpSHR = kron(identity(DR), Sp)
      SmSHR = kron(identity(DR), Sm)
      



      H = kron(HL, identity(DR * dim))  + kron(identity(DL * dim), HR) +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR)


      # Find Ground state



      (GSE,), GS = eigsh(H, k=1, which="SA")


      print (" ")
      print ("LN = ", LN )
      print ("E/L = ", GSE/LN)
      print (" ")




      # Construct reduced density matrix

     
      
      if HLMAT[p - 3].shape[0] < M:
          sd = HLMAT[p - 2].shape[0]  # p -2 because 16 x 16 matrix must still be computed before truncation
      else:
          sd = dim*M                    # 3 instead of 2 beacause the extra site is 3 dimensional not 2 dimensional.


     # print("HLMAT = ", HLMAT[p - 3].shape[0])

      psi = GS.reshape([sd, -1], order='C') 

    

      RDM = np.dot(psi, psi.conjugate().transpose())
      evalsl, evecsl = np.linalg.eigh(RDM)

     
      # Form the truncated reduced density matrix

      idx = evalsl.argsort()
      idx = idx[::-1] 



      if HLMAT[p - 2].shape[0] < M:
          m = HLMAT[p - 2].shape[0]
      else:
          m = M




      U = evecsl[:,idx[:m]]


      UENVL.append(U)         # Storing matrices for wave function transformation


     
    
      # Growth of operators for next iteration      

      if HRMAT[p  - 3].shape[0] < M:
          DR = HRMAT[p - 3].shape[0]
          DL = HRMAT[p - 3].shape[0]
      else:
          DR = m
          DL = m


      #left

      SztL = kron(identity(DL), Sz) 
      SptL = kron(identity(DL), Sp)
      SmtL = kron(identity(DL), Sm)
          
      # right 
         
      SztR = kron(identity(DR), Sz)
      SptR = kron(identity(DR), Sp) 
      SmtR = kron(identity(DR), Sm)






##############################################################################

# FINITE ALGORITHM



# Use HL from the last part of algorithm to grow left side.
# If you are ending with L = 20, it means that the system size is 10 sites. The following then will grow the left side by one site to make the system 11 sites. Therefore the environment is 9 sites. In turn this means that one should return to the case of L = 18 and store the right hand side block.   


# Start with symmetric configuration from infinite algorithm ---------00---------
# Grow left system by one site and shrink right system by one site. Keep in mind that the length of the chain remains at 20 sites.






M = [10, 20, 30, 40, 40]

ss = [30, 60, 90, 120, 120]   # Before, elements in M were multiplied by 2 to reflect that the additional site had 2-dimensional Hilbert; now it is 3-dimensional => now mulitply elements of M by 3 

# Refer to p.1714 for the following: 

# Bottom line: when the the loop indexed with j starts a new iteration the first run in the first loop indexed with i will have the same threshold of states, m, as the previous run. It is only AFTER the first run in the loop indexed with i that one truncates to a new m value. For example, in the above, m = 10 on the first run for j == 0. When one reaches the start of j == 1 run, m is still 10 at the beginning, but then m =20 is introduced towards of the first iteration in i for j == 1.  



# len( ) returns the length of the list.

for j in range (0, len(M) ):

    
 


# The following will end with ------------------00

   

    for i in range(1,int(It2) + 1):  # + 1 included because you are using "range", which starts at 1 and ends at int(It2)

        print("i = ", i)
       

        DR = HRMAT[-i-1].shape[0]/dim # dim = 3 instead of 2 since each site is 3 dimensional  



        if j != 0:

            DR = HRMAT[-i -  (int(It4) + 1)].shape[0]/dim   # Refer to p. 1763 as to why one has (int(It4) +  1)




      
      



        # Transformation of spin operators connecting left block to right block. 
        SztL = U.conjugate().transpose().dot(SztL.dot(U))  
        SptL = U.conjugate().transpose().dot(SptL.dot(U))
        SmtL = U.conjugate().transpose().dot(SmtL.dot(U))




        # Transformation of left system block

        HL = U.conjugate().transpose().dot(HL.dot(U)) 
      

        # Growth of left side by one site
           
        HL = kron(HL, identity(dim)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 

        HLMAT.append(HL)

        DL = HL.shape[0]/dim
     
       


         # left
  
        SzSHL = kron(identity(DL), Sz)
        SpSHL = kron(identity(DL), Sp)
        SmSHL = kron(identity(DL), Sm)

     


         # right
      
        SzSHR = kron(identity(DR), Sz)
        SpSHR = kron(identity(DR), Sp)
        SmSHR = kron(identity(DR), Sm)






# i == int(It2) signifies that the right system is now only composed of two sites. At this stage the right side becomes the system. In this current form of the code, this means that the superblock Hamilto      nian must be rearranged. Refer to p. 1700



        if i == int(It2):

           
           if j != 0:


               H = kron(HRMAT[-i -  (int(It4) + 1)], identity(DL * dim)) + kron(identity(DR * dim), HL)  +  0.5 * (kron(SpSHR, SmSHL)) +  0.5 * (kron(SmSHR, SpSHL)) + kron(SzSHR, SzSHL)

           elif j == 0:    

               H = kron(HRMAT[-i-1], identity(DL * dim)) + kron(identity(DR * dim), HL)  +  0.5 * (kron(SpSHR, SmSHL)) +  0.5 * (kron(SmSHR, SpSHL)) + kron(SzSHR, SzSHL)

        
        elif j != 0:

            H = kron(HL, identity(DR * dim)) + kron(identity(DL * dim), HRMAT[-i -  (int(It4)+ 1)])  +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR) 

        else:

            H = kron(HL, identity(DR * dim)) + kron(identity(DL * dim), HRMAT[-i-1])  +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR) 



       
   
        if j == 0:
            (GSE,), GS = eigsh(H, k=1, which="SA")
        else:
            (GSE,), GS = eigsh(H, k=1, which="SA", v0 = psi_trial)






        print (" ")
        print ("LN = ", LN )
        print ("E/L = ", GSE/Sit) # Sit is used because the lattice only possesses Sit sites.
        print (" ")




        

        if i == 2:
           sd = ss[j]


        if i == int(It2): # This will always be the case, no matter what value of sd is used
            sd = dim**2




        # Construct reduced density matrix


        psi = GS.reshape([sd, -1], order='C')

        RDM = np.dot(psi, psi.conjugate().transpose())
        evalsl, evecsl = np.linalg.eigh(RDM)


      
        # Contruct truncated transformation matrix

        idx = evalsl.argsort()
        idx = idx[::-1] 


   

        if i == int(It2):
                m = dim**2       # When dealing with spin - 1 system, this should be 9. In general, you could have dim**2 
        else:
                m = M[j]



  
        U = evecsl[:,idx[:m]]  






        if  i != int(It2):

            UENVL.append(U)         # Storing matrices for wave function transformation. != in place to take into account that the system and environment switch at this stage.

        elif j != 0 and i == int(It2):
            
            UENVR = []
            UENVR.append(U)
        
        else:

            UENVR.append(U)       # Because when i == int(It2) the left hand side becomes the environment and the right hand side becomes the system. With that, one stores the very last  





        ######## Break on final loop




        if j == 0 and  i == It2: # This is condition for wavefunction transformation to begin; it corresponds to configuration ------------------00


            psi_a = psi.reshape((-1, DL * dim), order="C")  # since the left hand side is the ennvironment dimension, DL is chosen here. AT THE MOMENT, I CANNOT EXPLAIN THIS PART
            psi_b = U.conjugate().transpose().dot(psi_a)
            psi_c = psi_b.reshape((-1, DL, dim), order="C").transpose(0, 2, 1)  # before, dim was 2. This works, but I do not know why.
            psi_d = psi_c.reshape((-1, DL), order="C")
            psi_trial = UENVL[-1].dot(psi_d.transpose()).transpose().reshape((-1, 1))   

           
           
           



        elif j != 0:   # Unlike the j == 0, where one  starts with ---------------00 configuration, j != 0 begins with ----------00----------, and the left hand side is the system until we reach the final iteratiion, when it swaps with the right hand side to become the environment. 


            if i == int (It2):


                psi_a = psi.reshape((-1, DL * dim), order="C") 
                psi_b = U.conjugate().transpose().dot(psi_a)
                psi_c = psi_b.reshape((-1, DL, dim), order="C").transpose(0, 2, 1)
                psi_d = psi_c.reshape((-1, DL), order="C")
                psi_trial = UENVL[-1].dot(psi_d.transpose()).transpose().reshape((-1, 1))   

            else:


                psi_a = psi.reshape((-1, DR * dim), order="C")  # since the left hand side is the ennvironment dimension, DL is chosen here. AT THE MOMENT, I CANNOT EXPLAIN THIS PART
                psi_b = U.conjugate().transpose().dot(psi_a)
                psi_c = psi_b.reshape((-1, DR, dim), order="C").transpose(0, 2, 1)
                psi_d = psi_c.reshape((-1, DR), order="C")
                psi_trial = UENVR[-i- (int(It4) + 1)].dot(psi_d.transpose()).transpose().reshape((-1, 1))    # -i - 9 because you wish to access the transformation associated with the configuration of one site extraon the left and one site less on the right. Also, the final element accessed in the last for loop of the previous is -9 => next element to be accessed is -10; since i begins with 1, -i - 9 leads immedaitely to -10. 



           
           
           




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

  
  




    for i in range(1, It3 + 1): # + 1 included to ensure that for loop goes through full range   

    
   

            DL = HLMAT[-i-1].shape[0]/dim # Refer to p. 1716; number of rows divided by 2 


            # Transformation of right side operators  

            SztR = U.conjugate().transpose().dot(SztR.dot(U)) 
            SptR = U.conjugate().transpose().dot(SptR.dot(U))
            SmtR = U.conjugate().transpose().dot(SmtR.dot(U))






            HR = U.conjugate().transpose().dot(HR.dot(U))   # HR instead of HL


             # Growth
           
            HR = kron(HR, identity(dim)) + kron(SztR, Sz) + 0.5 *( kron(SptR, Sm) + kron(SmtR, Sp)) 

            HRMAT.append(HR)
                
            DR = HRMAT[i].shape[0]/dim   # HRMAT list was empty just before this loop, with that the right hand system is grown and then stored in HRMAT using append. DR can then be obtained in a similar fashion                                 to the previous stage of the finite algorithm. However, instead of -i - 1, one has i since HRMAT is storing the matrices as the right hand system grows.
   


            # left
  
            SzSHL = kron(identity(DL), Sz) 
            SpSHL = kron(identity(DL), Sp)
            SmSHL = kron(identity(DL), Sm)

     

            # right
      

            SzSHR = kron(identity(DR), Sz) 
            SpSHR = kron(identity(DR), Sp)
            SmSHR = kron(identity(DR), Sm)






            if i == It3:

                H = kron(HLMAT[-i-1], identity(DR * dim)) + kron(identity(DL * dim), HR)  +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR)

            else:
  
                H = kron(HR, identity(DL * dim)) + kron(identity(DR * dim), HLMAT[-i-1])  +  0.5 * (kron(SpSHR, SmSHL)) +  0.5 * (kron(SmSHR, SpSHL)) + kron(SzSHR, SzSHL)






            # Find ground state       


            (GSE,), GS = eigsh(H, k=1, which="SA", v0 = psi_trial)




            print(" ")
            print(GSE/Sit)
            print(" ")
    





           # Construct reduced density matrix

       
            if HRMAT[i].shape[0] < dim* M[j]:   # it is 2*M instead of M beacuse you wish to compute 16x16 reduced density matrix before truncaiton. Now, dim instead of 2.
                        sd = HRMAT[i].shape[0]  
            elif i == int(It3):
                        sd = dim**2
            else:
                       sd = dim*M[j]
    



            psi = GS.reshape([sd, -1], order='C')    # Problem here when truncaiton begins. If statement needs to be implemented in order to ensure that the wavefunction is reshaped correctly.


            
            RDM = np.dot(psi, psi.conjugate().transpose())

            evalsl, evecsl = np.linalg.eigh(RDM)

            idx = evalsl.argsort()
            idx = idx[::-1] 







           # Form truncated matrix

            if HRMAT[i].shape[0] < M[j]: 
                        m = HRMAT[i].shape[0] 
            else:
                        m = M[j]



            # Right system truncation matrix. Needs to be stored in future, not now; still testing.             

            U = evecsl[:,idx[:m]] 

            if i != int(It3):

                UENVR.append(U)         # Storing matrices for wave function transformation. != in place to take into account that the system and environment switch at this stage.





            if i == int(It3):  # refer to p. 1761 of notes
                UENVL = []
                UENVL.append(U)



          






###########################

            if j == 0 and i < 7:    # This represents generating psi_trial up to approximately the middle of the chain. I am not sure why one does this in the simpe dmrg code. However, when one is on the second sweep, this apparently is no longer applicable. NOT SURE WHY THIS IS.
  
                psi_a = psi.reshape((-1, DL * dim), order="C")  # since the left hand side is the ennvironment dimension, DL is chosen here. AT THE MOMENT, I CANNOT EXPLAIN THIS PART
                psi_b = U.conjugate().transpose().dot(psi_a)
                psi_c = psi_b.reshape((-1, DL, dim), order="C").transpose(0, 2, 1)
                psi_d = psi_c.reshape((-1, DL), order="C")
                psi_trial = UENVL[-i - 1].dot(psi_d.transpose()).transpose().reshape((-1, 1))    # UENV[-i - 1] because we now wish to access -2 in the first iteration, -3 in the second iteration, etc.




            elif j == 0 and  i == int(It3):  # take into account that the system switches at this point in the iteration.
                     

                psi_a = psi.reshape((-1, DR * dim), order="C") 
                psi_b = U.conjugate().transpose().dot(psi_a)
                psi_c = psi_b.reshape((-1, DR, dim), order="C").transpose(0, 2, 1)
                psi_d = psi_c.reshape((-1, DR), order="C")      
                psi_trial = UENVR[-1].dot(psi_d.transpose()).transpose().reshape((-1, 1))  


            else:

                psi_trial = None







            if j != 0 and i != int(It3):   # Now we do not skip over any parts: no longer have psi_trial = None.

  
                psi_a = psi.reshape((-1, DL * dim), order="C")
                psi_b = U.conjugate().transpose().dot(psi_a)
                psi_c = psi_b.reshape((-1, DL, dim), order="C").transpose(0, 2, 1)
                psi_d = psi_c.reshape((-1, DL), order="C")
                psi_trial = UENVL[-i - 1].dot(psi_d.transpose()).transpose().reshape((-1, 1)) 



            elif j != 0 and i == int(It3):


                psi_a = psi.reshape((-1, DR * dim), order="C") 
                psi_b = U.conjugate().transpose().dot(psi_a)
                psi_c = psi_b.reshape((-1, DR, dim), order="C").transpose(0, 2, 1)
                psi_d = psi_c.reshape((-1, DR), order="C")
                psi_guess = UENVR[-1].dot(psi_d.transpose()).transpose().reshape((-1, 1))  




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







    
   


    for i in range(1, int(It4) + 1): # +1 included to allow for full range in loop   

            

                SztL = U.conjugate().transpose().dot(SztL.dot(U)) 
                SptL = U.conjugate().transpose().dot(SptL.dot(U))
                SmtL = U.conjugate().transpose().dot(SmtL.dot(U))



                HL = U.conjugate().transpose().dot(HL.dot(U))   


                # Growth
           
                HL = kron(HL, identity(dim)) + kron(SztL, Sz) + 0.5 *( kron(SptL, Sm) + kron(SmtL, Sp)) 

    
                HLMAT.append(HL)


                DL =  HLMAT[i].shape[0]/dim   # Refer to p. 1716; number of rows divided by 2        
    


                # left
  
                SzSHL = kron(identity(DL), Sz) 
                SpSHL = kron(identity(DL), Sp)
                SmSHL = kron(identity(DL), Sm)

     
                # right
    
                SzSHR = kron(identity(DR), Sz) 
                SpSHR = kron(identity(DR), Sp)
                SmSHR = kron(identity(DR), Sm)
   
    

                #Find HRMAT. HRNEW below is the HR matrix in the previous stage of this algorithm; it is the matrix associated with the penultimate run. 

                H = kron(HL, identity(DR * dim)) + kron(identity(DL * dim), HRMAT[-i-1]) +  0.5 * (kron(SmSHL, SpSHR)) +  0.5 * (kron(SpSHL, SmSHR)) + kron(SzSHL, SzSHR)




                # Find ground state

                (GSE,), GS = eigsh(H, k=1, which="SA", v0 = psi_trial)



                print(" ")
                print(GSE/Sit)
                print(" ")



                # Construct reduced density matrix


       
                if HLMAT[i].shape[0] < dim* M[j]:   # it is 2*M instead of M beacuse you wish to compute 16x16 reduced density matrix before truncaiton
                            sd = HLMAT[i].shape[0]  
                else:
                            sd = dim*M[j]




                psi = GS.reshape([sd, -1], order='C')    # Problem here when truncaiton begins. If statement needs to be implemented in order to ensure that the wavefunction is reshaped correctly.
                RDM = np.dot(psi, psi.conjugate().transpose())

                evalsl, evecsl = np.linalg.eigh(RDM)




                idx = evalsl.argsort()
                idx = idx[::-1] 




                # Form truncated matrix

                if HLMAT[i].shape[0] < M[j]: 
                            m = HLMAT[i].shape[0] 

                else:
                            m = M[j]


    


                U = evecsl[:,idx[:m]] 
                UENVL.append(U)






                # left
  
                SztL = kron(identity(DL), Sz) 
                SptL = kron(identity(DL), Sp)
                SmtL = kron(identity(DL), Sm)


                # right block. Refer to. 1676
             
                SztR = kron(identity(DR), Sz)
                SptR = kron(identity(DR), Sp) 
                SmtR = kron(identity(DR), Sm)

                

                if j == len(M) - 1 and i == int(It4): # in place as this ensures that the wave function transformation does not go one step too many.
                    break

                             

                psi_a = psi.reshape((-1, DR * dim), order="C") 
                psi_b = U.conjugate().transpose().dot(psi_a) 
                psi_c = psi_b.reshape((-1, DR, dim), order="C").transpose(0, 2, 1)         
                psi_d = psi_c.reshape((-1, DR), order="C")
                psi_trial = UENVR[-i -1].dot(psi_d.transpose()).transpose().reshape((-1, 1))    # UENV[-i - 1] because we now wish to access -2 in the first iteration, -3 in the second iteration, etc.




        






#####################################################
# The following creates the bond evolution operator for Suzuki - Trotter time evolution. 


H = [[0.25, 0, 0, 0], [0, -0.25, 0.5, 0], [0, 0.5 , -0.25, 0 ],[0, 0, 0, 0.25]]


GSE, GS = np.linalg.eigh(H)

####### INCLUDE FACTOR OF 0.5 FOR SECOND ORDER SUZUKI TROTTER

##### Make exponential matrix for both HA and HB. THERE IS A FACTOR OF 0.5 WITH HA





EV = GS.transpose().conjugate().dot(np.diag(np.exp(np.diag(GSE))).dot(GS))     # Time evolution operator for a single bond





print (time.time() - start_time, " seconds")
