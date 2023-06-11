*=======================================================================
*                                                                       
*       Author: Andreas Stathopoulos, Charlotte F. Fischer              
*                                                                       
*       Computer Science Department                                     
*       Vanderbilt University                                           
*       Nashville, TN 37212                                             
*       andreas@vuse.vanderbilt.edu                                     
*       cff@vuse.vanderbilt.edu                       DECEMBER 1993     
*                                                                       
*       Copyright (c) by Andreas Stathopoulos and Charlotte F. Fischer  
*                                                                       
*       DVDSON is a Fortran77 program that finds a few selected         
*       eigenvalues and their eigenvectors at either end of spectrum of 
*       a large, symmetric (and usually sparse) matrix, denoted as A.   
*       The matrix A is only referenced indirectly through the user     
*       supplied routine OP which implements a block matrix-vector      
*       operation(see below). Either the range of the eigenvalues wanted
*       or an array of the indices of selected ones can be specified.   
*       DVDSON is a front-end routine for setting up arrays, and initial
*       guess (calling SETUP). It also performs detailed error checking.
*       DVDRVR is the driver routine that implements a version of the   
*       Davidson algorithm. The characteristics of this version are:    
*        o  All arrays used by the program are stored in MEMORY.        
*        o  BLOCK method (many vectors may be targeted per iteration.)  
*        o  Eigenvectors are targeted in an optimum way without         
*           the need to compute all unconverged residuals,              
*        o  It REORTHOGONILIZES the basis in case of orthogonality loss.
*        o  Finds HIGHEST eigenpairs by using the negative of the A.    
*        o  Finds SELECTED eigenpairs specified by the user.            
*        o  It accepts INITIAL eigenvector ESTIMATES or it can          
*           CREATE INITIAL ESTIMATES from the diagonal elements.        
*        o  It uses a USER SUPPLIED block matrix-vector operation, OP.  
*           Depending on the implementation, OP can operate in either   
*           memory or on disc, and for either sparse or dense matrix.   
*        o  The user can provide STOPPING CRITERIA for eigenvalues,     
*           and residuals. The user can also CONTROL reorthogonalization
*            and block size.                                            
*        o  On exit INFORMATION is given about the convergence status   
*           of eigenpairs and the number of loops and OP operations.    
*                                                                       
*       The program consists of the following routines:                 
*       DVDSON, SETUP, DVDRVR, ADDABS, TSTSEL,                          
*       MULTBC, OVFLOW, NEWVEC, ORTHNRM.                                
                                                                        
*       It also calls some basic BLAS routines:                         
*       DCOPY, DSCAL, DDOT, DAXPY, IDAMAX, DGEMV, DINIT                 
                                                                        
*       For solving the small eigenproblem, the routine DSPEVX from     
*       LAPACK is used. DSPEVX is obtainable from NETLIB, together      
*       with a series of subroutines that it calls.                     
*                                                                       
*     All the routines have IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*-----------------------------------------------------------------------
*  (Important to the following is the concept of NUME, the distance of  
*   the index of the eigenpair wanted which is farthest from the        
*   extremes,i.e.,                                                      
*      if  lowest  eigepairs i1<i2<...<ik are wanted, NUME=ik           
*      if highest eigenpairs i1<i2<...<ik are wanted, NUME=N-i1+1       
*   where i1,...,ik are the indices of the wanted eigenpairs.           
*   Obviously, NUME.GE.(No. of EiGenpairs wanted). )                    
                                                                        
*   on entry                                                            
*   -------                                                             
*   OP          User supplied routine with calling sequence OP(N,M,B,C).
*               B and C are N x M matrices and C stores the result AxB. 
*               It should be declared external in the main program.     
*   N           Order of the matrix.                                    
*   LIM         The upper limit on the dimension of the expanding basis.
*               NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed
*               only for LIM=NUME=N. The choice of LIM depends on the   
*               available workspace (see below). If the space is        
*               available it is preferable to have a large LIM, but not 
*               larger than NUME$+$40.                                  
*   DIAG        Array of size N with the diagonal elements of the       
*               matrix A.                                               
*   ILOW        The index of the lowest eigepair to be computed. If     
*               (ILOW.LE.0).or.(ILOW.GT.N), the selected eigenpairs     
*               to be computed should be contained in array ISELEC.     
*               (Modified on exit).                                     
*   IHIGH       The index of the highest eigenpair to be computed.      
*               Considered ONLY when ILOW is in the range               
*               (0.LT.ILOW.LE.N). (Modified on exit).                   
*   ISELEC      Array of size LIM holding the user specified indices    
*               for the eigenpairs to be computed. Considered only when 
*               (ILOW.LE.0).or.(ILOW.GT.N). The indices are read from   
*               the first position until a non positive integer is met. 
*                  Example: if N=500, ILOW=0, and ISELEC(1)=495,        
*                  ISELEC(2)=497, ISELEC(3)=-1, the program will find   
*                  2 of the highest eigenpairs, pairs 495 and 497.      
*               Any order of indices is acceptable (Modified on exit).  
*   NIV         Number of Initial Vector estimates provided by the user.
*               If NIV is in the range:  (NUME).LE.(NIV).LE.(LIM),      
*               the first NIV columns of size N of WORK should contain  
*               the estimates (see below). In all other cases of NIV,   
*               the program generates initial estimates.                
*   MBLOCK      Number of vectors to be targeted in each iteration.     
*               1.LE.MBLOCK.LE.(No. EiGenpairs wanted) should hold.     
*               Large block size reduces the number of iterations       
*               (matrix acceses) but increases the matrix-vector        
*               multiplies. It should be used when the matrix accese    
*               is expensive (disc, recomputed or distributed).         
*   CRITE       Convergence threshold for eigenvalues.                  
*               If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted 
*               eigenvalues, convergence is signaled.                   
*   CRITC       Convergence threshold for the coefficients of the last  
*               added basis vector(s). If all of those corresponding to 
*               unconverged eigenpairs are less than CRITC convergence  
*               is signaled.                                            
*   CRITR       Convergence threshold for residual vector norms. If     
*               all the residual norms ||Ax_i-l_ix_i|| of the targeted  
*               x_i are less than CRITR convergence is signaled.        
*               If ANY of the criteria are satisfied the algorithm stops
*   ORTHO       The threshold over which loss of orthogonality is       
*               assumed. Usually ORTHO.LE.CRITR*10 but the process can  
*               be skipped by setting ORTHO to a large number(eg,1.D+3).
*   MAXITER     Upper bound on the number of iterations of the          
*               algorithm. When MAXITER is exceeded the algorithm stops.
*               A typical MAXITER can be MAX(200,NUME*40), but it can   
*               be increased as needed.                                 
*   WORK        Real array of size IWRSZ. Used for both input and output
*               If NIV is in ((NUME).LE.(NIV).LE.(LIM)), on input, WORK 
*               must have the NIV initial estimates. These NIV N-element
*               vectors start from WORK(1) and continue one after the   
*               other. They must form an orthonormal basis.             
*   IWRSZ       The size of the real workspace. It must be at least as  
*               large as:                                               
*                                                                       
*                       2*N*LIM + LIM*LIM + (NUME+10)*LIM + NUME        
*                                                                       
*   IWORK       Integer work array of size IIWSZ. Used as scrath array  
*               for indices and for use in the LAPACK routines.         
*   IIWSZ       The size of the integer workspace. It must be at least  
*               as large as:                                            
*                                    6*LIM + NUME                       
*                                                                       
*               If LIM or NUME needs to be increased, the space should  
*               also be increased accordingly. For given IWRSZ and      
*               IIWSZ one can calculate how big a problem one can       
*               solve (LIM,NUME).                                       
*                                                                       
*   on exit                                                             
*   -------                                                             
*   WORK(1)     The first NUME*N locations contain the approximations to
*               the NUME extreme eigenvectors. If the lowest eigenpairs 
*               are required, (HIEND=false), eigenvectors appear in     
*               ascending order, otherwise (HIEND=false), they appear in
*               descending order. If only some are requested, the order 
*               is the above one for all the NUME extreme eigenvectors, 
*               but convergence has been reached only for the selected  
*               ones. The rest are the current approximations to the    
*               non-selected eigenvectors.                              
*   WORK(NUME*N+1)                                                      
*               The next NUME locations contain the approximations to   
*               the NUME extreme eigenvalues, corresponding to the above
*               NUME eigenvectors. The same ordering and convergence    
*               status applies here as well.                            
*   WORK(NUME*N+NUME+1)                                                 
*               The next NUME locations contain the corresponding values
*               of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of 
*               the last step of the algorithm.                         
*   WORK(NUME*N+NUME+NUME+1)                                            
*               The next NUME locations contain the corresponding       
*               residual norms of the NUME above eigenvectors, of the   
*               last step.                                              
*   HIEND       Logical. If .true. on exit the highest eigenpairs are   
*               found in descending order. Otherwise, the lowest        
*               eigenpairs are arranged in ascending order.             
*   NLOOPS      The number of iterations it took to reach convergence.  
*               This is also the number of matrix references.           
*   NMV         The number of Matrix-vector(M-V) multiplies. Each matrix
*               reference can have up to size(block) M-V multiplies.    
*   IERR        An integer denoting the completions status:             
*               IERR = 0        denotes normal completion.              
*               IERR = -k       denotes error in DSPEVX (k eigenpairs   
*                               not converged)                          
*               0<IERR<=2048    denotes some inconsistency as follows:  
*        If (INT( MOD(IERR,  2)/1  ) N < LIM                            
*        If (INT( MOD(IERR,  4)/2  ) LIM < 1                            
*        If (INT( MOD(IERR,  8)/4  ) ISELEC(1)<1, and no range specified
*        If (INT( MOD(IERR, 16)/8  ) IHIGH > N (in range or ISELEC)     
*        If (INT( MOD(IERR, 32)/16 ) IHIGH < ILOW (Invalid range)       
*        If (INT( MOD(IERR, 64)/32 ) NEIG >= LIM (Too many wanted)      
*        If (INT( MOD(IERR,128)/64 ) Probable duplication in ISELEC     
*        If (INT( MOD(IERR,256)/128) NUME >= LIM (max eigen very far)   
*        If (INT( MOD(IERR,512)/256) MBLOCK is out of bounds            
*        If (INT( MOD(IERR,1024)/512) IWRSZ or IIWSZ is not enough      
*        If (INT( MOD(IERR,2048)/1024) Orthogonalization Failed         
*        If (INT( MOD(IERR,4096)/2048) NLOOPS > MAXITER                 
*                                                                       
*               The program will also print an informative message to   
*               the standard output when NIV is not proper but it will  
*               continue by picking initial estimates internally.       
      SUBROUTINE DVDSON(OP,N,LIM,DIAG,                                
     :     ILOW,IHIGH,ISELEC,NIV,MBLOCK,                  
     :     CRITE,CRITC,CRITR,ORTHO,MAXITER,               
     :     WORK,IWRSZ,IWORK,IIWSZ,                        
     :     HIEND,NLOOPS,NMV,IERR)                         
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
      DIMENSION DIAG(N),WORK(IWRSZ),IWORK(IIWSZ)                      
      DIMENSION ISELEC(LIM)                                           
      LOGICAL HIEND                                                   
      EXTERNAL OP                                                     
      INTEGER IERR
*-----------------------------------------------------------------------
*     
*     Checking user input errors, and setting up the problem to solve.      
*     
      IERR=0                                                          
      IF (LIM.GT.N) IERR=IERR+1                                       
      IF (LIM.LE.0) IERR=IERR+2                                       
      
      HIEND=.false.                                                   
      
      IF ((ILOW.LE.0).OR.(ILOW.GT.N)) THEN                            
*     ..Look for user choice of eigenpairs in ISELEC               
         IF (ISELEC(1).LE.0) THEN                                     
*     ..Nothing is given in ISELEC                              
            IERR=IERR+4                                               
         ELSE                                                         
*     ..Find number of eigenpairs wanted, and their             
*     ..min/max indices                                         
            NEIG=1                                                    
            ILOW=ISELEC(1)                                            
            IHIGH=ISELEC(1)                                           
            DO I=2,LIM                                             
               IF (ISELEC(I).LE.0) GOTO 20                            
               ILOW=MIN(ILOW,ISELEC(I))                               
               IHIGH=MAX(IHIGH,ISELEC(I))                             
               NEIG=NEIG+1                                            
            ENDDO
*     ..Check if a very large index is asked for                
 20         IF (IHIGH.GT.N) IERR=IERR+8                               
         ENDIF                                                        
      ELSE                                                            
*     ..Look for a range between ILOW and IHIGH                    
*     ..Invalid range. IHIGH>N                                     
         IF (IHIGH.GT.N) IERR=IERR+8                                  
         NEIG=IHIGH-ILOW+1                                            
*     ..Invalid range. IHIGH<ILOW                                  
         IF (NEIG.LE.0) IERR=IERR+16                                  
         IF (NEIG.GT.LIM) THEN                                        
*     ..Not enough Basis space. Increase LIM or decrease NEIG   
            IERR=IERR+32                                              
         ELSE                                                         
*     ..Fill in the ISELEC with the required indices            
            DO I=1,NEIG                                               
               ISELEC(I)=ILOW+I-1                                     
            ENDDO
         ENDIF                                                        
      ENDIF                                                           
      
      IF (IERR.NE.0) RETURN                                           
      
      NUME=IHIGH                                                      
*     ..Identify if few of the highest eigenpairs are wanted.         
      IF ((ILOW+IHIGH-1).GT.N) THEN                                   
         HIEND=.true.                                                 
         NUME=N-ILOW+1                                                
*     ..Change the problem to a minimum eipenpairs one             
*     ..by picking the corresponding eigenpairs on the             
*     ..opposite side of the spectrum.                             
         DO I=1,NEIG                                                  
            ISELEC(I)=N-ISELEC(I)+1
         ENDDO
      ENDIF                                                           
*     ..duplications in ISELEC                                        
      IF (NEIG.GT.NUME) IERR=IERR+64                                  
*     ..Not enough Basis space. Increase LIM or decrease NUME         
      IF ((NUME.GT.LIM).OR.((NUME.EQ.LIM).AND.(NUME.NE.N)))           
     :     IERR=IERR+128                                                
*     .     .Size of Block out of bounds                                   
      IF ( (MBLOCK.LT.1).OR.(MBLOCK.GT.NEIG) ) IERR=IERR+256          
      
*     ..Check for enough workspace for Dvdson                         
      IF ((IWRSZ.LT.(LIM*(2*N+LIM+(NUME+10))+NUME)).OR.               
     :     (IIWSZ.LT.(6*LIM+NUME))) IERR=IERR+512                      
      
      IF (IERR.NE.0) RETURN                                           
      
      IF (NIV.GT.LIM) THEN                                            
*     ..Check number of initial estimates NIV is lower than LIM.   
         PRINT*,'WARNING: Too many initial estimates.?'               
         PRINT*,'The routine will pick the appropriate number'        
      ELSEIF ((NIV.LT.NUME).AND.(NIV.GT.0)) THEN                      
*     ..check if enough initial estimates.                         
*     ..(NIV<1 => program chooses)                                 
         PRINT*,'WARNING: Not enough initial estimates'               
         PRINT*,'The routine will pick the appropriate number'        
      ENDIF                                                           
*     
*     Assigning space for the real work arrays                              
*     
      iBasis    =1                                                    
      ieigval   =iBasis  +N*LIM                                       
      iAB       =ieigval +LIM                                         
      iS        =iAB     +N*LIM                                       
      itempS    =iS      +LIM*(LIM+1)/2                               
      iSvec     =itempS  +LIM*(LIM+1)/2                               
      iscra1    =iSvec   +LIM*NUME                                    
      ioldval   =iscra1  +8*LIM                                       
*     
*     Assigning space for the integer work arrays                           
*     
      iscra2    =1                                                    
      iscra3    =iscra2  +5*LIM                                       
      iIcv      =iscra3  +LIM                                         
      
      IF (HIEND) CALL DSCAL(N,-1.D0,DIAG,1)                           
      
      iSTART=NIV                                                      
      CALL SETUP(OP,N,LIM,NUME,HIEND,DIAG,WORK(iscra1),               
     :     WORK(iBasis),WORK(iAB),WORK(iS),iSTART)
      NLOOPS=1                                                        
      NMV=ISTART                                                      
      
      CALL DVDRVR(OP,N,HIEND,LIM,MBLOCK,DIAG,                         
     :     NUME,iSTART,NEIG,ISELEC,                             
     :     CRITE,CRITC,CRITR,ORTHO,MAXITER,                     
     :     WORK(ieigval),WORK(iBasis),WORK(iAB),                
     :     WORK(iS),WORK(itempS),WORK(iSvec),                   
     :     WORK(iscra1),IWORK(iscra2),IWORK(iscra3),            
     :     IWORK(iIcv),WORK(ioldval),                           
     :     NLOOPS,NMV,IERR)                                     
      
      IF (HIEND) THEN                                                 
         CALL DSCAL(N,-1.D0,DIAG,1)                                   
         CALL DSCAL(NUME,-1.D0,WORK(ieigval),1)                       
      endif                                                           
*     
*     -Copy the eigenvalues after the eigenvectors                          
*     -Next, copy the difference of eigenvalues between the last two steps  
*     -Next, copy the residuals for the first NUME estimates                
*     
      CALL DCOPY(NUME,WORK(ieigval),1,WORK(iBasis+N*NUME),1)          
      CALL DCOPY(NUME,WORK(ioldval),1,WORK(iBasis+(N+1)*NUME),1)      
      CALL DCOPY(NUME,WORK(iscra1),1,WORK(iBasis+(N+2)*NUME),1)       
      
 100  RETURN                                                          
      END                                                             
*=======================================================================
      SUBROUTINE SETUP(OP,N,LIM,NUME,HIEND,DIAG,MINELEM,              
     :     BASIS,AB,S,NIV)                                
*=======================================================================
*     Subroutine for setting up (i) the initial BASIS if not provided,
*     (ii) the product of the matrix A with the Basis into matrix AB, 
*     and (iii) the small matrix S=B^TAB. If no initial estimates are 
*     available, the BASIS =(e_i1,e_i2,...,e_iNUME), where i1,i2,..., 
*     iNUME are the indices of the NUME lowest diagonal elements, and 
*     e_i the i-th unit vector. (ii) and (iii) are handled by ADDABS. 
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
      DIMENSION DIAG(N),BASIS(N*LIM),AB(N*LIM)                        
      DIMENSION S(LIM*(LIM+1)/2),MINELEM(LIM)
      EXTERNAL OP
      LOGICAL HIEND
      DOUBLE PRECISION MINELEM
*-----------------------------------------------------------------------
*     on entry                                                            
*     --------                                                            
*     OP          The block matrix-vector operation, passed to ADDABS     
*     N           the order of the matrix A                               
*     LIM         The limit on the size of the expanding Basis            
*     NUME        Largest index of the wanted eigenvalues.                
*     HIEND       Logical. True only if the highest eigenpairs are needed.
*     DIAG        Array of size N with the diagonal elements of A         
*     MINELEM     Array keeping the indices of the NUME lowest diagonals. 
*     
*     on exit                                                             
*     -------                                                             
*     BASIS       The starting basis.                                     
*     AB, S       The starting D=AB, and small matrix S=B^TAB             
*     NIV         The starting dimension of BASIS.                        
*-----------------------------------------------------------------------
      
      IF ((NIV.GT.LIM).OR.(NIV.LT.NUME)) THEN                         
*     
*     ..Initial estimates are not available. Give as estimates unit
*     ..vectors corresponding to the NUME minimum diagonal elements
*     ..First find the indices of these NUME elements (in MINELEM).
*     ..Array AB is used temporarily as a scratch array.           
*     
         CALL DINIT(N,-1.D0,AB,1)                                     
         DO I=1,NUME                                               
*     ..imin= the first not gotten elem( NUME<=N )              
            DO 20 J=1,N                                               
 20            IF (AB(J).LT.0) GOTO 30                                
 30            IMIN=J                                                    
               DO 40 J=IMIN+1,N                                          
 40               IF ((AB(J).LT.0).AND.                                  
     :                 (DIAG(J).LT.DIAG(IMIN))) IMIN=J                     
                  MINELEM(I)=IMIN                                           
                  AB(IMIN)=1.D0                                             
               ENDDO
*     
*     ..Build the Basis. B_i=e_(MINELEM(i))                        
*     
               CALL DINIT(N*LIM,0.D0,BASIS,1)                               
               DO J=1,NUME                                               
                  I=(J-1)*N+MINELEM(J)                                      
                  BASIS(I)=1                                                
               ENDDO
               
               NIV=NUME                                                     
            ENDIF                                                           
*     
*     Find the matrix AB by matrix-vector multiplies, as well as the        
*     small matrix S = B^TAB.                                               
*     
            KPASS=0                                                         
            CALL ADDABS(OP,N,LIM,HIEND,KPASS,NIV,BASIS,AB,S)                
            
            RETURN                                                          
            END                                                             
*=======================================================================
      SUBROUTINE DVDRVR(OP,N,HIEND,LIM,MBLOCK,DIAG,                   
     :     NUME,NIV,NEIG,ISELEC,                         
     :     CRITE,CRITC,CRITR,ORTHO,MAXITER,              
     :     EIGVAL,BASIS,AB,S,TEMPS,SVEC,                 
     :     SCRA1,ISCRA2,INCV,ICV,OLDVAL,                 
     :     NLOOPS,NMV,IERR)                              
*=======================================================================
*     called by DVDSON                                                
*     
*     Driver routine implementing Davidson's main loop. On entry it   
*     is given the Basis, the work matrix D=AB and the small symmetric
*     matrix to be solved, S=B^TAB (as found by SETUP). In each step  
*     the small problem is solved by calling DSPEVX.                  
*     TSTSEL tests for eigenvalue convergence and selects the next    
*     pairs to be considered for targeting (as a block).              
*     NEWVEC computes the new vectors (block) to be added in the      
*     expanding basis, and tests for residual convergence.            
*     ADDABS is the critical step of matrix multiplication. The new   
*     vectors of D are found Dnew=ABnew, and the new small problem S, 
*     is calculated. The algorithm is repeated.                       
*     In case of a large expanding basis (KPASS=LIM) the Basis, AB,   
*     SVEC and S are collapsed.                                       
*     At the end the current eigenvector estimates are computed as    
*     well as the residuals and eigenvalue differences.               
*     
*     Subroutines called:                                             
*     DSPEVX, MULTBC, TSTSEL, OVFLOW, NEWVEC, ADDABS,                 
*     DCOPY, DDOT, DAXPY                                              
*-----------------------------------------------------------------------
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
      DIMENSION DIAG(N)                                               
      DIMENSION S(LIM*(LIM+1)/2),TEMPS(LIM*(LIM+1)/2)                 
      DIMENSION SVEC(LIM*NUME),EIGVAL(LIM)                            
      DIMENSION ISELEC(NEIG)                                          
      DIMENSION BASIS(N*LIM),AB(N*LIM)                                
      DIMENSION SCRA1(8*LIM),ISCRA2(5*LIM),INCV(LIM)                  
      DIMENSION ICV(NUME),OLDVAL(NUME)                                
      LOGICAL RESTART,FIRST,DONE,HIEND,TSTSEL                         
      EXTERNAL OP                                                     
*-----------------------------------------------------------------------
*     
*     on entry                                                            
*     -------                                                             
*     
*     OP          The user specified block-matrix-vector routine          
*     N           The order of the matrix A                               
*     HIEND       Logical. True only if the highest eigenpairs are needed.
*     LIM         The limit on the size of the expanding Basis            
*     MBLOCK      Number of vectors to be targeted in each iteration.     
*     DIAG        Array of size N with the diagonal elements of A         
*     NUME        The largest index of the eigenvalues wanted.            
*     NIV         Starting dimension of expanding basis.                  
*     NEIG        Number of eigenvalues wanted.                           
*     ISELEC      Array containg the indices of those NEIG eigenpairs.    
*     CRITE       Convergence thresholds for eigenvalues, coefficients    
*     CRITC,CRITR and residuals.                                          
*     BASIS       Array with the basis vectors.                           
*     AB          Array with the vectors D=AB                             
*     S           Array keeping the symmetric matrix of the small problem.
*     TEMPS       scratch array                                           
*     SVEC        Array for holding the eigenvectors of S                 
*     SCRA1       Srcatch array used by DSPEVX.                           
*     ISCRA2      Integer Srcatch array used by DSPEVX.                   
*     INCV        Srcatch array used in DSPEVX. Also used in TSTSEL and   
*     NEWVEC where it holds the Indices of uNConVerged pairs  
*     ICV         It contains "1" to the locations of ConVerged eigenpairs
*     OLDVAL      Array keeping the previous' step eigenvalue estimates.  
*     
*     on exit                                                             
*     -------                                                             
*     
*     EIGVAL      Array containing the NUME lowest eigenvalues of the     
*     the matrix A (or -A if the highest are sought).         
*     Basis       On exit Basis stores the NUME corresponding eigenvectors
*     OLDVAL      On exit it stores the final differences of eigenvalues. 
*     SCRA1       On exit it stores the NUME corresponding residuals.     
*     NLOOPS      Number of loops taken by the algorithm                  
*     NMV         Number of matrix-vector products performed.             
*     
*-----------------------------------------------------------------------
      DO I=1,NUME                                                   
         EIGVAL(I)=1.D30                                              
         ICV(I)=0
      ENDDO
      FIRST =.true.                                                   
      KPASS =NIV                                                      
      NNCV  =KPASS                                                    
      
 10   CONTINUE                                                        
*     (iterations for kpass=NUME,LIM)                                 
*     
*     Diagonalize the matrix S. Find only the NUME smallest eigenpairs      
*     
      CALL DCOPY(NUME,EIGVAL,1,OLDVAL,1)                           
      CALL DCOPY((KPASS*(KPASS+1))/2,S,1,TEMPS,1)                  
      CALL DSPEVX('Vectors also','In a range','Upper triangular',  
     :     KPASS,TEMPS,-1.,-1.,1,NUME,0.D0,                        
     :     NFOUND,EIGVAL,SVEC,KPASS,SCRA1,ISCRA2,INCV,INFO)        
      IERR=-ABS(INFO)                                              
      IF (IERR.NE.0) GOTO 60                                       
*     
*     TeST for convergence on the absolute difference of eigenvalues between
*     successive steps. Also SELect the unconverged eigenpairs and sort them
*     by the largest magnitude in the last added NNCV rows of Svec.         
*     
      DONE=TSTSEL(KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV,          
     :     CRITE,CRITC,SCRA1,ISCRA2,OLDVAL,NNCV,INCV)            
      IF ((DONE).OR.(KPASS.GE.N)) GOTO 30                          
      
      IF (KPASS.EQ.LIM) THEN                                       
*     Maximum size for expanding basis. Collapse basis, D, and S, Svec      
*     Consider the basis vectors found in TSTSEL for the newvec.            
*     
         CALL MULTBC(N,LIM,NUME,SVEC,SCRA1,BASIS)                  
         CALL MULTBC(N,LIM,NUME,SVEC,SCRA1,AB)                     
         CALL OVFLOW(NUME,LIM,S,SVEC,EIGVAL)                       
         KPASS=NUME                                                
      ENDIF                                                        
*     
*     Compute and add the new vectors. NNCV is set to the number of new     
*     vectors that have not converged. If none, DONE=true, exit.            
*     
      CALL NEWVEC(N,NUME,LIM,MBLOCK,KPASS,CRITR,ORTHO,NNCV,INCV,   
     :     DIAG,SVEC,EIGVAL,AB,BASIS,ICV,RESTART,DONE)      
      
*     ..An infinite loop is avoided since after a collapsing Svec=I
*     ..=> Res=Di-lBi which is just computed and it is orthogonal. 
*     ..The following is to prevent an improbable infinite loop.   
      IF (.NOT.RESTART) THEN                                       
         FIRST=.true.                                              
      ELSEIF (FIRST) THEN                                          
         FIRST=.false.                                             
         CALL MULTBC(N,KPASS+NNCV,NUME,SVEC,SCRA1,BASIS)           
         CALL MULTBC(N,KPASS+NNCV,NUME,SVEC,SCRA1,AB)              
         CALL OVFLOW(NUME,KPASS+NNCV,S,SVEC,EIGVAL)                
         KPASS=NUME                                                
         GOTO 10                                                   
      ELSE                                                         
         IERR=IERR+1024                                            
         GOTO 30                                                   
      ENDIF                                                        
      
      IF (DONE) GOTO 30                                            
*     
*     Add new columns in D and S, from the NNCV new vectors.                
*     
      CALL ADDABS(OP,N,LIM,HIEND,KPASS,NNCV,BASIS,AB,S)            
      
      NMV=NMV+NNCV                                                 
      KPASS=KPASS+NNCV                                             
      NLOOPS=NLOOPS+1                                              
      
      IF (NLOOPS.LE.MAXITER) GOTO 10                                  
      IERR=IERR+2048                                                  
      NLOOPS=NLOOPS-1                                                 
      KPASS=KPASS-NNCV                                                
 30   CONTINUE                                                        
*     
*     Calculate final results. EIGVAL contains the eigenvalues, BASIS the   
*     eigenvectors, OLDVAL the eigenvalue differences, and SCRA1 residuals. 
*     
      DO I=1,NUME                                                  
         OLDVAL(I)=ABS(OLDVAL(I)-EIGVAL(I))
      ENDDO
      
      CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,BASIS)                      
      CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,AB)                         
*     
*     i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi                        
*     temporarily stored in AB(NUME*N+1)                                    
*     
      DO I=1,NUME                                                  
         CALL DCOPY(N,AB((I-1)*N+1),1,AB(NUME*N+1),1)                 
         CALL DAXPY(N,-EIGVAL(I),BASIS((I-1)*N+1),1,AB(NUME*N+1),1)   
         SCRA1(I)=DDOT(N,AB(NUME*N+1),1,AB(NUME*N+1),1)               
         SCRA1(I)=SQRT(SCRA1(I))                                      
      ENDDO
      
 60   RETURN                                                          
      END                                                             
*=======================================================================
        SUBROUTINE ADDABS(OP,N,LIM,HIEND,KPASS,NNCV,BASIS,AB,S)         
*=======================================================================
*       Called by: DVDRVR, SETUP                                        
*                                                                       
*       Calculates the new column in the D matrix and the new column    
*       in the S matrix. The new D column is D(new)=AB(new). S has a    
*       new row and column, but being symmetric only the new column is  
*       stored. S(i,kpass+1)=B(i)^T D(kpass+1) for all i.               
*                                                                       
*       subroutines called:                                             
*       OP, DDOT, DSCAL                                                 
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
        DIMENSION BASIS(N*LIM),AB(N*LIM)                                
        DIMENSION S(LIM*(LIM+1)/2)                                      
        LOGICAL HIEND                                                   
        EXTERNAL OP                                                     
*-----------------------------------------------------------------------
*   on entry                                                            
*   -------                                                             
*   N           The order of the matrix A                               
*   kpass       The current dimension of the expanding sub-basis        
*   NNCV        Number of new basis vectors.                            
*   Basis       the basis vectors, including the new NNCV ones.         
*   on exit                                                             
*   -------                                                             
*   AB          The new matrix D=AB. (with new NNCV columns)            
*   S           The small matrix with NNCV new columns at the last part 
*-----------------------------------------------------------------------
*                                                                       
* The user specified matrix-vector routine is called with the new       
* basis vector B(*,kpass+1) and the result is assigned to AB(idstart)   
*                                                                       
        IDSTART=KPASS*N+1                                               
        CALL OP(N,NNCV,BASIS(IDSTART),AB(IDSTART))                      
*                                                                       
* If highest pairs are sought, use the negative of the matrix           
*                                                                       
        IF (HIEND) CALL DSCAL(N*NNCV,-1.D0,AB(IDSTART),1)               
*                                                                       
* The new S is calculated by adding the new last columns                
* S(new)=B^T D(new).                                                    
*                                                                       
        ISSTART=KPASS*(KPASS+1)/2                                       
        DO 20 IV=1,NNCV                                                 
           IBSTART=1                                                    
           DO 10 IBV=1,KPASS+IV                                         
               SS=DDOT(N,BASIS(IBSTART),1,AB(IDSTART),1)                
               S(ISSTART + IBV)=SS                                      
               IBSTART=IBSTART+N                                        
 10        CONTINUE                                                     
           ISSTART=ISSTART+KPASS+IV                                     
           IDSTART=IDSTART+N                                            
 20     CONTINUE                                                        
                                                                        
        RETURN                                                          
        END                                                             
*=======================================================================
        LOGICAL FUNCTION TSTSEL(KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV, 
     :                       CRITE,CRITC,ROWLAST,IND,OLDVAL,NNCV,INCV)  
*=======================================================================
*                                                                       
*       Called by: DVDRVR                                               
                                                                        
*       It first checks if the wanted eigenvalues have reached          
*       convergence and updates OLDVAL. Second, for each wanted and non 
*       converged eigenvector, it finds the largest absolute coefficient
*       of the NNCV last added vectors (from SVEC) and if not coverged, 
*       places it in ROWLAST. IND has the corresponding indices.        
*       Third, it sorts ROWLAST in decreasing order and places the      
*       corresponding indices in the array INCV. The index array INCV   
*       and the number of unconverged pairs NNCV, are passed to DVDRVR. 
*       Later in NEWVEC only the first MBLOCK of NNCV pairs will be     
*       targeted, since if ROWLAST(i) > ROWLAST(j)                      
*       then approximately RESIDUAL(i) > RESIDUAL(j)                    
*                                                                       
*       Subroutines called                                              
*       IDAMAX                                                          
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
        LOGICAL DONE                                                    
        DIMENSION SVEC(KPASS*NUME),EIGVAL(NUME)                         
        DIMENSION ICV(NUME)                                             
        DIMENSION ROWLAST(NEIG),IND(NEIG),OLDVAL(NUME)                  
        DIMENSION INCV(NEIG),ISELEC(NEIG)                               
*-----------------------------------------------------------------------
*                                                                       
*   on entry                                                            
*   -------                                                             
*   KPASS       current dimension of the expanding Basis                
*   NUME        Largest index of the wanted eigenvalues.                
*   NEIG        number of wanted eigenvalues of original matrix         
*   ISELEC      index array of the wanted eigenvalues.                  
*   SVEC        the eigenvectors of the small system                    
*   EIGVAL      The NUME lowest eigenvalues of the small problem        
*   ICV         Index of converged eigenpairs.ICV(i)=1 iff eigenpair i  
*               has converged, and ICV(i)=0 if eigenpair i has not.     
*   CRITE,CRITC Convergence thresholds for eigenvalues and coefficients 
*   ROWLAST     scratch array, keeping the largest absolute coefficient 
*               of the NNCV last rows of Svec.                          
*   IND         scratch array, temporary keeping the indices of Rowlast 
*   OLDVAL      The previous iteration's eigenvalues.                   
*                                                                       
*   on exit                                                             
*   -------                                                             
*   NNCV         Number of non converged eigenvectors (to be targeted)  
*   INCV         Index to these columns in decreasing order of magnitude
*   TSTSEL       true if convergence has been reached                   
*                                                                       
*-----------------------------------------------------------------------
                                                                        
        DONE=.False.                                                    
*                                                                       
* Test all wanted eigenvalues for convergence under CRITE               
*                                                                       
        NNCE=0                                                          
        DO 10 I=1,NEIG                                                  
           IVAL=ISELEC(I)                                               
 10        IF (ABS(OLDVAL(IVAL)-EIGVAL(IVAL)).GE.CRITE) NNCE=NNCE+1     
        IF (NNCE.EQ.0) THEN                                             
           TSTSEL=.TRUE.                                                
           RETURN                                                       
        ENDIF                                                           
*                                                                       
* Find the maximum element of the last NNCV coefficients of unconverged 
* eigenvectors. For those unconverged coefficients, put their indices   
* to IND and find their number NNCV                                     
*                                                                       
        ICNT=0                                                          
        DO I=1,NEIG                                                  
           IF (ICV(ISELEC(I)).EQ.0) THEN                                
*             ..Find coefficient and test for convergence               
              ICUR=KPASS*ISELEC(I)                                      
              TMAX=ABS( SVEC(ICUR) )                                    
              DO 20 L=1,NNCV-1                                          
 20              TMAX=MAX( TMAX, ABS(SVEC(ICUR-L)) )                    
              IF (TMAX.LT.CRITC) THEN                                   
*                ..this  coefficient converged                          
                 ICV(ISELEC(I))=1                                       
              ELSE                                                      
*                ..Not converged. Add it to the list.                   
                 ICNT=ICNT+1                                            
                 IND(ICNT)=ISELEC(I)                                    
                 ROWLAST(ICNT)=TMAX                                     
              ENDIF                                                     
           ENDIF                                                        
        ENDDO
        
        NNCV=ICNT                                                       
        IF (NNCV.EQ.0) DONE=.TRUE.                                      
*                                                                       
* Sort the ROWLAST elements interchanging their indices as well         
*                                                                       
        DO 40 I=1,NNCV                                                  
           INDX=IDAMAX(NNCV-I+1,ROWLAST(I),1)                           
           INCV(I)=IND(INDX+I-1)                                        
                                                                        
           TEMP=ROWLAST(INDX+I-1)                                       
           ROWLAST(INDX+I-1)=ROWLAST(I)                                 
           ROWLAST(I)=TEMP                                              
           ITEMP=IND(INDX+I-1)                                          
           IND(INDX+I-1)=IND(I)                                         
           IND(I)=ITEMP                                                 
 40     CONTINUE                                                        
                                                                        
        TSTSEL=DONE                                                     
        RETURN                                                          
        END                                                             
*=======================================================================
        SUBROUTINE MULTBC(N,K,M,C,TEMP,B)                               
*=======================================================================
*       called by: DVDRVR                                               
*                                                                       
*       Multiplies B(N,K)*C(K,M) and stores it in B(N,M)                
*       Used for collapsing the expanding basis to current estimates,   
*       when basis becomes too large, or for returning the results back 
                                                                        
*       Subroutines called                                              
*       DINIT, DGEMV, DCOPY                                             
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
        DIMENSION B(N*K),C(K*M),TEMP(M)                                 
*-----------------------------------------------------------------------
        DO 10 IROW=1,N                                                  
*              CALL DINIT(M,0.d0,TEMP,1)                                
           CALL DGEMV('Transp',K,M, 1.D0, C,K,B(IROW),N, 0.D0 ,TEMP,1)  
           CALL DCOPY(M,TEMP,1,B(IROW),N)                               
  10    CONTINUE                                                        
                                                                        
        RETURN                                                          
        END                                                             
                                                                        
*=======================================================================
        SUBROUTINE OVFLOW(NUME,LIM,S,SVEC,EIGVAL)                       
*=======================================================================
*       Called by: DVDRVR                                               
*       Called when the upper limit (LIM) has been reached for the basis
*       expansion. The new S is computed as S'(i,j)=l(i)delta(i,j) where
*       l(i) eigenvalues, and delta of Kronecker, i,j=1,NUME. The new   
*       eigenvectors of the small matrix are the unit vectors.          
*                                                                       
*       Subroutines called:                                             
*       DCOPY, DINIT                                                    
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
        DIMENSION SVEC(LIM*NUME),S((LIM*(LIM+1))/2)                     
        DIMENSION EIGVAL(LIM)                                           
*-----------------------------------------------------------------------
*   on entry                                                            
*   -------                                                             
*   NUME        The largest index of eigenvalues wanted.                
*   SVEC        the kpass eigenvectors of the smaller system solved     
*   EIGVAL      the eigenvalues of this small system                    
*   on exit                                                             
*   -------                                                             
*   S           The new small matrix to be solved.                      
*-----------------------------------------------------------------------
*                                                                       
* calculation of the new upper S=diag(l1,...,l_NUME) and                
* its matrix Svec of eigenvectors (e1,...,e_NUME)                       
*                                                                       
        CALL DINIT((NUME*(NUME+1))/2,0.d0,S,1)                          
        CALL DINIT(NUME*NUME,0.d0,SVEC,1)                               
        IND=0                                                           
        ICUR=0                                                          
        DO 10 I=1,NUME                                                  
           S(IND+I)=EIGVAL(I)                                           
           SVEC(ICUR+I)=1                                               
           ICUR=ICUR+NUME                                               
   10      IND=IND+I                                                    
                                                                        
        RETURN                                                          
        END                                                             
*=======================================================================
        SUBROUTINE NEWVEC(N,NUME,LIM,MBLOCK,KPASS,CRITR,ORTHO,NNCV,INCV,
     :                    DIAG,SVEC,EIGVAL,AB,BASIS,ICV,RESTART,DONE)   
*=======================================================================
*                                                                       
*       Called by: DVDRVR                                               
*                                                                       
*       It calculates the new expansion vectors of the basis.           
*       For each one of the vectors in INCV starting with the largest   
*       megnitude one, calculate its residual Ri= DCi-liBCi and check   
*       the ||Ri|| for convergence. If it is converged do not add it    
*       but look for the immediate larger coefficient and its vector.   
*       The above procedure continues until MBLOCK vectors have been    
*       added to the basis, or the upper limit has been encountered.    
*       Thus only  the required MBLOCK residuals are computed. Then,    
*       calculate the first order correction on the added residuals     
*       Ri(j) = Ri(j)/(li-Ajj) and orthonormalizes the new vectors      
*       to the basis and to themselves.                                 
*                                                                       
*       Subroutines called:                                             
*       ORTHNRM, DDOT, DGEMV                                            
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
        DIMENSION INCV(NUME)                                            
        DIMENSION ICV(NUME)                                             
        DIMENSION DIAG(N)                                               
        DIMENSION BASIS(N*LIM),AB(N*LIM)                                
        DIMENSION SVEC(LIM*NUME)                                        
        DIMENSION EIGVAL(LIM)                                           
        LOGICAL RESTART,DONE                                            
*-----------------------------------------------------------------------
*   on entry                                                            
*   --------                                                            
*   N           The order of the matrix A                               
*   NUME        The largest index of the eigenvalues wanted.            
*   LIM         The limit on the size of the expanding Basis            
*   MBLOCK      Maximum number of vectora to enter the basis            
*   KPASS       the current dimension of the expanding basis            
*   CRITR       Convergence threshold for residuals                     
*   ORTHO       Orthogonality threshold to be passed to ORTHNRM         
*   NNCV        Number of Non ConVerged pairs (MBLOCK will be targeted) 
*   INCV        Index to the corresponding SVEC columns of these pairs. 
*   DIAG        Array of size N with the diagonal elements of A         
*   SVEC,EIGVAL Arrays holding the eigenvectors and eigenvalues of S    
*   AB          Array with the vectors D=AB                             
*   BASIS       the expanding basis having kpass vectors                
*   ICV         Index of converged eigenpairs (ICV(i)=1 <=>i converged) 
                                                                        
*   on exit                                                             
*   -------                                                             
*   NNCV        The number of vectors finally added to the basis.       
*   BASIS       The new basis incorporating the new vectors at the end  
*   ICV         Index of converged eigenpairs (updated)                 
*   DONE        logical, if covergance has been reached.                
*   RESTART     logical, if because of extreme loss of orthogonality    
*               the Basis should be collapsed to current approximations.
*-----------------------------------------------------------------------
        DONE    = .FALSE.                                               
        NEWSTART= KPASS*N+1                                             
        NADDED  = 0                                                     
        ICVC    = 0                                                     
        LIMADD  = MIN( LIM, MBLOCK+KPASS )                              
        ICUR    = NEWSTART                                              
*                                                                       
* Compute RESIDUALS for the MBLOCK of the NNCV not converged vectors.   
*                                                                       
        DO 10 I=1,NNCV                                                  
           INDX=INCV(I)                                                 
*          ..Compute  Newv=BASIS*Svec_indx , then                       
*          ..Compute  Newv=AB*Svec_indx - eigval*Newv and then          
*          ..compute the norm of the residual of Newv                   
           CALL DGEMV('N',N,KPASS,1.D0,BASIS,N,SVEC((INDX-1)*KPASS+1),1,
     :                 0.d0,BASIS(ICUR),1)                              
           CALL DGEMV('N',N,KPASS,1.D0,AB,N,SVEC((INDX-1)*KPASS+1),1,   
     :                 -EIGVAL(INDX),BASIS(ICUR),1)                     
           SS = DNRM2(N,BASIS(ICUR),1)                                  
*                                                                       
*          ..Check for convergence of this residual                     
*                                                                       
           IF (SS.LT.CRITR) THEN                                        
*             ..Converged,do not add. Go for next non converged one     
              ICVC=ICVC+1                                               
              ICV( INDX ) = 1                                           
              IF (ICVC.LT.NNCV) GOTO 10                                 
*             ..All have converged.                                     
              DONE=.TRUE.                                               
              RETURN                                                    
           ELSE                                                         
*             ..Not converged. Add it in the basis                      
              NADDED=NADDED+1                                           
              INCV(NADDED)=INDX                                         
              IF ((NADDED+KPASS).EQ.LIMADD) GOTO 20                     
*             ..More to be added in the block                           
              ICUR=ICUR+N                                               
           ENDIF                                                        
 10     CONTINUE                                                        
                                                                        
 20     NNCV=NADDED                                                     
*                                                                       
* Diagonal preconditioning: newvect(i)=newvect(i)/(l-Aii)               
* If (l-Aii) is very small (or zero) divide by 10.D-6                   
*                                                                       
        ICUR=NEWSTART-1                                                 
        DO 50 I=1,NNCV                                                  
           DO 40 IROW=1,N                                               
              DG=EIGVAL(INCV(I))-DIAG(IROW)                             
              IF (ABS(DG).GT.(1.D-13)) THEN                             
                  BASIS(ICUR+IROW)=BASIS(ICUR+IROW) / DG                
              ELSE                                                      
                  BASIS(ICUR+IROW)=BASIS(ICUR+IROW) /1.D-13             
              ENDIF                                                     
 40        CONTINUE                                                     
           ICUR=ICUR+N                                                  
 50     CONTINUE                                                        
*                                                                       
* ORTHONORMALIZATION                                                    
*                                                                       
        CALL ORTHNRM(N,LIM,ORTHO,KPASS,NNCV,AB(NEWSTART),               
     :               BASIS,RESTART)                                     
                                                                        
 99     RETURN                                                          
        END                                                             
*=======================================================================
        SUBROUTINE ORTHNRM(N,LIM,ORTHO,KPASS,NNCV,SCRA1,                
     :                     BASIS,RESTART)                               
*=======================================================================
*                                                                       
*       It orthogonalizes the new NNCV basis vectors starting from the  
*       kpass+1, to the previous vectors of the basis and to themselves.
*       A Gram-Schmidt method is followed after which the residuals     
*       should be orthogonal to the BASIS. Because of machine arithmetic
*       errors this orthogonality may be lost, and a reorthogonalization
*       procedure is adopted whenever orthogonality loss is above a     
*       ORTHO. If after some reorthogonalizations the procedure does not
*       converge to orthogonality, the basis is collapsed to the        
*       current eigenvector approximations.                             
*                                                                       
*       Subroutines called:                                             
*       DAXPY, DDOT, DSCAL                                              
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
        DIMENSION BASIS(N*LIM)                                          
        DIMENSION SCRA1(N)                                              
        LOGICAL RESTART                                                 
*-----------------------------------------------------------------------
*   on entry                                                            
*   --------                                                            
*   N           The order of the matrix A                               
*   LIM         The limit on the size of the expanding Basis            
*   ORTHO       The orthogonality threshold                             
*   KPASS       The number of basis vectors already in Basis            
*   NNCV        The number of new vectors in the basis                  
*   SCRA1       Scratch vector of size N                                
*   BASIS       the expanding basis having kpass vectors                
*                                                                       
*   on exit                                                             
*   -------                                                             
*   BASIS       the new basis orthonormalized                           
*   RESTART     Logical, if true the algoritm will collapse BASIS.      
*-----------------------------------------------------------------------
*                                                                       
* ORTHOGONALIZATION                                                     
*                                                                       
        RESTART=.false.                                                 
        ICUR=KPASS*N+1                                                  
*                                                                       
*       .. do iv=1,nncv                                                 
        IV = 1                                                          
 30     CONTINUE                                                        
                                                                        
           DPREV=1.D+7                                                  
 5         DCUR=0.D0                                                    
           IBSTART=1                                                    
           DO 10 I=1,KPASS+IV-1                                         
              SCRA1(I)=DDOT(N,BASIS(IBSTART),1,BASIS(ICUR),1)           
              DCUR=MAX(DCUR,ABS(SCRA1(I)))                              
              IBSTART=IBSTART+N                                         
 10        CONTINUE                                                     
           IBSTART=1                                                    
           DO 20 I=1,KPASS+IV-1                                         
              CALL DAXPY(N,-SCRA1(I),BASIS(IBSTART),1,BASIS(ICUR),1)    
              IBSTART=IBSTART+N                                         
 20        CONTINUE                                                     
                                                                        
           IF (DCUR.GE.ORTHO) THEN                                      
              IF (DCUR.GT.DPREV) THEN                                   
                 RESTART=.true.                                         
*                ..Adjust the number of added vectors.                  
                 NNCV=IV-1                                              
                 RETURN                                                 
              ELSE                                                      
                 DPREV=DCUR                                             
                 GOTO 5                                                 
              ENDIF                                                     
           ENDIF                                                        
*                                                                       
* NORMALIZATION                                                         
*                                                                       
           SCRA1(1)=DDOT(N,BASIS(ICUR),1,BASIS(ICUR),1)                 
           SCRA1(1)=SQRT(SCRA1(1))                                      
           IF (SCRA1(1).LT.1D-14) THEN                                  
              CALL DCOPY(N,BASIS( N*(NNCV-1)+1),1,BASIS(ICUR),1)        
              NNCV=NNCV-1                                               
           ELSE                                                         
              CALL DSCAL(N,1/SCRA1(1),BASIS(ICUR),1)                    
              ICUR=ICUR+N                                               
              IV = IV +1                                                
           ENDIF                                                        
        IF (IV.LE.NNCV) GOTO 30                                         
                                                                        
        RETURN                                                          
        END
      
c$$$        PROGRAM DRIVER_DAV                                              
c$$$************************************************************************
c$$$* This driver program is to demonstrate several ways on how the routine 
c$$$* DVDSON should be called, what input it expects and what output it     
c$$$* gives. This program by no means exhausts all the different input      
c$$$* and output combinations of the algorithm, but it gives four useful    
c$$$* examples.                                                             
c$$$* A lower triangual matrix is created and for the first two calls the   
c$$$* 10 lowest eigenvalues are computed. These two examples demonstrate    
c$$$* how to specify a range of eigenvalues for computation, how to vary    
c$$$* the convergence criteria, and how to give initial estimates by setting
c$$$* NIV and the first NIV vectors of WORK.                                
c$$$* The next two exmamples compute highest eigenpairs. The first one      
c$$$* demonstrates how to selectively choose only a few eigenpairs to be    
c$$$* computed, while both of them show how the upper extreme is handled.   
c$$$* The importance of reorthogonilization to some cases is also shown.    
c$$$*                                                                       
c$$$* Important in this program is the notion of NUME (see DVDSON comments).
c$$$* NUME is the index of the farthest from the extreme computed eigenvalue
c$$$* DVDSON always comnute approximations from the extreme to this NUME,   
c$$$* even if it was asked for the last only. However these approximations  
c$$$* will not be accurate except the ones explicitly selected.             
c$$$*                                                                       
c$$$* Typical values for the thresholds and controls are:                   
c$$$*                                                                       
c$$$*  CRITE < from 1D-10 to 1D-14                                          
c$$$*  CRITC < from 1D-9 to 1D-14                                           
c$$$*  CRITR < from 1D-6 to 1D-10                                           
c$$$*  Of course these can be varied according to the specific needs.       
c$$$*  Eigenvalues converge much faster than residuals, so a very accurate  
c$$$*  computation would require CRITE <1D-15 and CRITR<1D-8. However       
c$$$*  residuals depend on the norm of A. Eigenvectors converging under     
c$$$*  |coef|<CRITC have components accurate to about log10(CRITC) digits.  
c$$$                                                                        
c$$$*  ORTHO < CRITR. (=1D+30 is off, =CRITR*10 is tightly on).             
c$$$*  If the user is confident about  the convergence characteristics      
c$$$*  of the problem reorthogonalization can be switched off by setting    
c$$$*  ORTHO >1D+10 (or something big). However, specific cases like the    
c$$$*  fourth of the cases below may need it when tight convergence is      
c$$$*  required (CRITR < 1D-8)                                              
c$$$                                                                        
c$$$*  MBLOCK (1<=MBLOCK<=No. of EiG. wanted)                               
c$$$*  How many of the unconverged vectors are targeted in each iteration.  
c$$$*  Block methods (large MBLOCK) ususally work well with matrix-vector   
c$$$*  multiplies for disc or for parallel machines or when the matrix is   
c$$$*  not stored but computed at each iteration. If the matrix is stored in
c$$$*  memory then block method can be switched off by setting MBLOCK=1.    
c$$$*                                                                       
c$$$*  MAXITER                                                              
c$$$*  It has been observed that DVDSON takes from 10-30 iterations         
c$$$*  per eigenpair. However setting MAXITER a little higher does not      
c$$$*  require any extra storage and it could save from rerunning the code  
c$$$*  for a few more iterations.                                           
c$$$                                                                        
c$$$************************************************************************
c$$$                                                                        
c$$$        IMPLICIT DOUBLE PRECISION (A-H,O-Z)                             
c$$$        PARAMETER(                                                      
c$$$     :            Nmax     =  100,                                      
c$$$     :            IBAND    =  10,                                       
c$$$     :            NZERmax  =  (IBAND+1)*(Nmax-IBAND/2),                 
c$$$     :            NUMEmax  =  10,                                       
c$$$     :            LIMmax   =  NUMEmax+20,                               
c$$$     :            IRWSZ    =  2*Nmax*LIMmax + LIMmax*LIMmax +           
c$$$     :                        (NUMEmax+10)*LIMmax + NUMEmax ,           
c$$$     :            IIWSZ    =  6*LIMmax + NUMEmax )                      
c$$$                                                                        
c$$$        COMMON /MATRIX/ A(NZERmax),IndCol(Nmax),IndRow(NZERmax),LUPPER  
c$$$        COMMON /TEMP/ TEMPB(Nmax), TEMPC(Nmax)                          
c$$$        DIMENSION DIAG(Nmax)                                            
c$$$        DIMENSION WORK(IRWSZ),IWORK(IIWSZ)                              
c$$$        DIMENSION ISELEC(LIMmax)                                        
c$$$        LOGICAL LUPPER                                                  
c$$$        EXTERNAL DSSBMV                                                 
c$$$*                                                                       
c$$$* Create a test-matrix(lower triangular) in the sparse format described 
c$$$* in the paper. It is banded with band symmetric band size = iband. All 
c$$$* off diagonals are set to 0.001 and the diagonals to their column index
c$$$*                                                                       
c$$$        N=Nmax                                                          
c$$$        ICUR=1                                                          
c$$$        DO 10 I=1,N                                                     
c$$$           A(ICUR)=I                                                    
c$$$           INDROW(ICUR)=I                                               
c$$$           INDCOL(I)=ICUR                                               
c$$$           ICUR=ICUR+1                                                  
c$$$           DO 20 J=I+1,MIN(N,I+IBAND)                                   
c$$$              A(ICUR)=0.001                                             
c$$$              INDROW(ICUR)=J                                            
c$$$              INDCOL(I)=ICUR                                            
c$$$ 20           ICUR=ICUR+1                                               
c$$$ 10     CONTINUE                                                        
c$$$        LUPPER=.False.                                                  
c$$$                                                                        
c$$$        DIAG(1)=A(1)                                                    
c$$$        DO 30 I=2,N                                                     
c$$$ 30        DIAG(I)=A( INDCOL(I-1)+1 )                                   
c$$$                                                                        
c$$$************************************************************************
c$$$*               TEsting                                                 
c$$$************************************************************************
c$$$*                                                                       
c$$$* Case 1. The lowest range [l(1),...,l(NUMEmax)] is asked, with no      
c$$$* initial estimates available. The threshold for residual is really big 
c$$$* so that approximate estimates can be found in a few steps.            
c$$$*                                                                       
c$$$        N       = Nmax                                                  
c$$$        LIM     = LIMmax                                                
c$$$        ILOW    = 1                                                     
c$$$        IHIGH   = NUMEmax                                               
c$$$        NUME    = NUMEmax                                               
c$$$        NIV     = 0                                                     
c$$$        CRITE   = 1D-15                                                 
c$$$        CRITC   = 1D-8                                                  
c$$$        CRITR   = 1D-3                                                  
c$$$        MBLOCK  = 1                                                     
c$$$        ORTHO   = 1D+3                                                  
c$$$        MAXITER = MAX( NUME*40, 200 )                                   
c$$$                                                                        
c$$$        CALL DVDSON(DSSBMV,N,LIM,DIAG,                                  
c$$$     :             ILOW,IHIGH,ISELEC,NIV,MBLOCK,                        
c$$$     :             CRITE,CRITC,CRITR,ORTHO,MAXITER,                     
c$$$     :             WORK,IRWSZ,IWORK,IIWSZ,HIEND,NLOOPS,NMV,IERR)        
c$$$                                                                        
c$$$*                                                                       
c$$$* OUTPUT                                                                
c$$$        WRITE(6,1000) IERR,NLOOPS,NMV,HIEND                             
c$$$ 1000   FORMAT(/'IERR =',I5,'  Matrix accesses=',I4,                    
c$$$     :   '  Matrix-Vector products=',I4,// '  Descending Order: ',L1)   
c$$$        WRITE(6,2000) ((WORK(i), i=NUME*N+j,(N+3)*NUME,NUME),j=1,NUME)  
c$$$ 2000   FORMAT(//9X,'Eigenvalues',8X,'Eigval Differences',6X,           
c$$$     :    'Residuals',//(D25.15,2D20.10))                               
c$$$        WRITE(6,3000) ((WORK(I), I=J,J+5), J=1,NUME*N,N)                
c$$$ 3000   FORMAT(//' First six componenents of the eigenvectors'//        
c$$$     :    (6D12.3))                                                     
c$$$                                                                        
c$$$************************************************************************
c$$$* Run again the same problem but now with tight residual threshold      
c$$$* and initial estimates the ones obtained above:                        
c$$$* The rest of the input variables remain the same.                      
c$$$*                                                                       
c$$$        NIV     = NUME                                                  
c$$$        CRITC   = 1D-12                                                 
c$$$        CRITR   = 1D-8                                                  
c$$$                                                                        
c$$$        CALL DVDSON(DSSBMV,N,LIM,DIAG,                                  
c$$$     :             ILOW,IHIGH,ISELEC,NIV,MBLOCK,                        
c$$$     :             CRITE,CRITC,CRITR,ORTHO,MAXITER,                     
c$$$     :             WORK,IRWSZ,IWORK,IIWSZ,HIEND,NLOOPS,NMV,IERR)        
c$$$                                                                        
c$$$*                                                                       
c$$$* OUTPUT                                                                
c$$$        WRITE(6,1000) IERR,NLOOPS,NMV,HIEND                             
c$$$        WRITE(6,2000) ((WORK(i), i=NUME*N+j,(N+3)*NUME,NUME),j=1,NUME)  
c$$$        WRITE(6,3000) ((WORK(I), I=J,J+5), J=1,NUME*N,N)                
c$$$                                                                        
c$$$************************************************************************
c$$$* CASE 2.                                                               
c$$$* Find only selected eigenvalues of the highest part of the spectrum:   
c$$$* Here we ask for N, N-5, N-9. No initial estimates are specified.      
c$$$* We allow for block method by setting MBLOCK=3 (no. of pairs wanted)   
c$$$* We also switch the reorthogonalization procedure over a threshold     
c$$$                                                                        
c$$$        ILOW    = 0                                                     
c$$$        ISELEC(1)= N-9                                                  
c$$$        ISELEC(2)= N                                                    
c$$$        ISELEC(3)= N-5                                                  
c$$$        ISELEC(4)= -1                                                   
c$$$* Notice the value of NUME=10. This is because the NUmber of Maximum    
c$$$* Eigenpair needed which is farthest from the highest extreme, is 10    
c$$$* for the N-9. DVDSON will give rough approximations to the rest of the 
c$$$* eigenpairs up to the N-9 (ten of them)                                
c$$$        NUME    = 10                                                    
c$$$        NIV     = 0                                                     
c$$$        CRITE   = 1D-14                                                 
c$$$        CRITC   = 1D-12                                                 
c$$$        CRITR   = 1D-7                                                  
c$$$        MBLOCK  = 3                                                     
c$$$        ORTHO   = 1D-6                                                  
c$$$        MAXITER = MAX( NUME*40, 200 )                                   
c$$$                                                                        
c$$$        CALL DVDSON(DSSBMV,N,LIM,DIAG,                                  
c$$$     :             ILOW,IHIGH,ISELEC,NIV,MBLOCK,                        
c$$$     :             CRITE,CRITC,CRITR,ORTHO,MAXITER,                     
c$$$     :             WORK,IRWSZ,IWORK,IIWSZ,HIEND,NLOOPS,NMV,IERR)        
c$$$******                                                                  
c$$$* Notice the decreasing order:                                          
c$$$******                                                                  
c$$$*                                                                       
c$$$* OUTPUT                                                                
c$$$        WRITE(6,1000) IERR,NLOOPS,NMV,HIEND                             
c$$$        WRITE(6,2000) ((WORK(i), i=NUME*N+j,(N+3)*NUME,NUME),j=1,NUME)  
c$$$        WRITE(6,3000) ((WORK(I), I=J,J+5), J=1,NUME*N,N)                
c$$$                                                                        
c$$$************************************************************************
c$$$* Use the above estimates as approximations to find all the eigenpairs  
c$$$* from N-9 to N (highest extreme). If reorthogonalization is switched   
c$$$* off, the current example cannot converge to the accuracy specified    
c$$$* below.                                                                
c$$$*                                                                       
c$$$                                                                        
c$$$        ILOW = N-9                                                      
c$$$        IHIGH= N                                                        
c$$$        NUME = 10                                                       
c$$$        NIV  = 10                                                       
c$$$        MBLOCK=10                                                       
c$$$        CRITE= 1D-14                                                    
c$$$        CRITC= 1D-12                                                    
c$$$        CRITR= 1D-10                                                    
c$$$        ORTHO= 1D-9                                                     
c$$$                                                                        
c$$$        CALL DVDSON(DSSBMV,N,LIM,DIAG,                                  
c$$$     :             ILOW,IHIGH,ISELEC,NIV,MBLOCK,                        
c$$$     :             CRITE,CRITC,CRITR,ORTHO,MAXITER,                     
c$$$     :             WORK,IRWSZ,IWORK,IIWSZ,HIEND,NLOOPS,NMV,IERR)        
c$$$******                                                                  
c$$$* Notice the decreasing order:                                          
c$$$******                                                                  
c$$$*                                                                       
c$$$* OUTPUT                                                                
c$$$        WRITE(6,1000) IERR,NLOOPS,NMV,HIEND                             
c$$$        WRITE(6,2000) ((WORK(i), i=NUME*N+j,(N+3)*NUME,NUME),j=1,NUME)  
c$$$        WRITE(6,3000) ((WORK(I), I=J,J+5), J=1,NUME*N,N)                
c$$$                                                                        
c$$$        STOP                                                            
c$$$        END                                                             
c$$$                                                                        
c$$$*=======================================================================
c$$$        SUBROUTINE DSSBMV(N,M,B,C)                                      
c$$$*=======================================================================
c$$$*                                                                       
c$$$*       Computes the product of matrix A with a block of vectors B(N,M) 
c$$$*                               C=A B                                   
c$$$*       where A(NxN) is a Symmetric Sparse matrix. Only the nonzero     
c$$$*       elements of either the upper or lower triangular matrix are     
c$$$*       kept and managed by Indices.                                    
c$$$*                                                                       
c$$$*       Subroutines called:                                             
c$$$*       DGATHR,DSCATR,DDOT,DAXPY,DINIT                                  
c$$$*=======================================================================
c$$$*                                                                       
c$$$        IMPLICIT DOUBLE PRECISION(A-H,O-Z)                              
c$$$        PARAMETER(Nmax=100,NZERmax=1045)                                
c$$$        COMMON/MATRIX/  A(NZERmax),INDCOL(Nmax),INDROW(NZERmax),LUPPER  
c$$$        COMMON/TEMP/ TEMPB(Nmax),TEMPC(Nmax)                            
c$$$        DIMENSION B(N*M),C(N*M)                                         
c$$$        LOGICAL LUPPER                                                  
c$$$                                                                        
c$$$************************************************************************
c$$$*                                                                       
c$$$*   on entry                                                            
c$$$*   --------                                                            
c$$$*   N           the order of the matrix A                               
c$$$*   B           the matrix (block of vectors) to multiply A with        
c$$$*   A           Linear array keeping the the nonzero elements of        
c$$$*               the matrix A. It stores columns one after the other     
c$$$*               starting from the first. It is either the upper         
c$$$*               triangular part or the lower part depending on logical  
c$$$*               LUPPER. (max elements that may contain= n^2/(2*nodes))  
c$$$*   INDCOL      It is the index showing for each column where that      
c$$$*               column ends in array A.                                 
c$$$*   INDROW      It is the index showing for each element of A to its    
c$$$*               row number in the square matrix                         
c$$$*   LUPPER      logical. If .true. the upper part of A is used.         
c$$$*   TEMPB,TEMPC Linear scratch arrays of size N                         
c$$$*                                                                       
c$$$*   on exit                                                             
c$$$*   -------                                                             
c$$$*                                                                       
c$$$*   C           the result of the multiplication (dim=NxM)              
c$$$*                                                                       
c$$$************************************************************************
c$$$        IOFFS=1                                                         
c$$$        IF (LUPPER) IOFFS=0                                             
c$$$                                                                        
c$$$        ISTART=1                                                        
c$$$        CALL DINIT(N*M,0.D0,C,1)                                        
c$$$                                                                        
c$$$        DO 20 ICOL=1,N                                                  
c$$$           NUMELEM = INDCOL(ICOL)-ISTART+1                              
c$$$           ICUR=1                                                       
c$$$                                                                        
c$$$           DO 10 IV=1,M                                                 
c$$$              CALL GATHER(NUMELEM,TEMPB,B(ICUR),INDROW(ISTART))         
c$$$              CALL GATHER(NUMELEM,TEMPC,C(ICUR),INDROW(ISTART))         
c$$$              DIAG=C(ICUR-1+ICOL)+DDOT(NUMELEM,A(ISTART),1,TEMPB,1)     
c$$$              CALL DAXPY(NUMELEM-1,B(ICUR-1+ICOL),                      
c$$$     :                    A(ISTART+IOFFS),1,                            
c$$$     :                    TEMPC(IOFFS+1),1)                             
c$$$                                                                        
c$$$              CALL SCATTER(NUMELEM,C(ICUR),INDROW(ISTART),TEMPC)        
c$$$              C(ICUR-1+ICOL)=DIAG                                       
c$$$              ICUR=ICUR+N                                               
c$$$   10      CONTINUE                                                     
c$$$                                                                        
c$$$            ISTART=INDCOL(ICOL)+1                                       
c$$$   20   CONTINUE                                                        
c$$$                                                                        
c$$$        RETURN                                                          
c$$$        END                                                             
c$$$                                                                        
c$$$*=======================================================================
c$$$        SUBROUTINE GATHER(N,A,B,INDEX)                                  
c$$$*=======================================================================
c$$$*       This subroutine collects array elements accessed via an         
c$$$*       integer pointer to contiguous storage.                          
c$$$*=======================================================================
c$$$                                                                        
c$$$        INTEGER N,INDEX(1)                                              
c$$$        DOUBLE PRECISION A(1),B(1)                                      
c$$$                                                                        
c$$$        DO 10 I = 1,N                                                   
c$$$          A(I) = B(INDEX(I))                                            
c$$$   10   CONTINUE                                                        
c$$$                                                                        
c$$$        RETURN                                                          
c$$$        END                                                             
c$$$                                                                        
c$$$*=======================================================================
c$$$        SUBROUTINE SCATTER(N,A,INDEX,B)                                 
c$$$*=======================================================================
c$$$*       This subroutine disperses array elements accessed via an integer
c$$$*       pointer from contiguous storage to the appropriate location.    
c$$$*=======================================================================
c$$$                                                                        
c$$$        INTEGER N,INDEX(1)                                              
c$$$        DOUBLE PRECISION A(1),B(1)                                      
c$$$                                                                        
c$$$        DO 10 I = 1,N                                                   
c$$$          A(INDEX(I)) = B(I)                                            
c$$$   10   CONTINUE                                                        
c$$$                                                                        
c$$$        RETURN                                                          
c$$$        END                                                             
                                                                        
*=======================================================================
        SUBROUTINE DINIT( N, A, X, INCX )                               
*=======================================================================
*       PURPOSE ... INITIALIZES DOUBLE PRECISION VECTOR TO              
*                   A CONSTANT VALUE 'A'                                
*=======================================================================
        INTEGER N, INCX                                                 
        DOUBLE PRECISION A, X (*)                                       
        INTEGER XADDR, I                                                
                                                                        
        IF  ( INCX .EQ. 1 )  THEN                                       
            DO 100 I = 1, N                                             
                X(I) = A                                                
  100       CONTINUE                                                    
                                                                        
        ELSE                                                            
            XADDR = 1                                                   
            IF  ( INCX .LT. 0 )  THEN                                   
                XADDR = (-N+1)*INCX + 1                                 
            ENDIF                                                       
            DO 200 I = 1, N                                             
                X(XADDR) = A                                            
                XADDR = XADDR + INCX                                    
  200       CONTINUE                                                    
        ENDIF                                                           
                                                                        
        RETURN                                                          
        END                                                             
                                                                        
