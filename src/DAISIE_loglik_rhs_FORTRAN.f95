
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! Example how to apply Fortran code with variable-length parameters
! ---------------------------------------------------------------------
! ------------------------------------------------------------------------

!==========================================================================
! Helper function: 
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE daisie_fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
      II = II
        DO I = 1, DIMP
          II = II + 1
          vec(I) = parms(II)
        ENDDO
        
      END SUBROUTINE daisie_fill1d

!==========================================================================
! module with declarations
!==========================================================================

      MODULE daisie_dimmod

      ! length of the vector -  decided in R-code
      INTEGER  :: N
      INTEGER  :: kk
      
      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)  
      
      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE daisie_dimmod

!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE daisie_initmod (steadyparms)
      USE daisie_dimmod 

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 2  ! constant-length parameters
      
      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block 

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector       
      N = INT(parms(1) + 1e-6)
      kk = INT(parms(2) + 1e-6)

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)  
      ALLOCATE(P(5 * (N + 4 + 2 * kk)))

      initialised = .FALSE.
       
      END SUBROUTINE daisie_initmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE daisie_runmod (neq, t, Conc, dConc, yout, ip)
      USE daisie_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(2 * N + 1), dConc(2 * N + 1), yout(*)
      DOUBLE PRECISION  :: xx1(N + 3), xx2(N + 3), xx3
      INTEGER 		      :: il1(N), il2(N), il3in3(N), il4(N)
	    INTEGER           :: in1(N), in2ix2(N)
      INTEGER           :: ix1(N), ix3(N), ix4(N)
      DOUBLE PRECISION  :: laavec(N + 4 + 2 * kk),lacvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: muvec(N + 4 + 2 * kk),gamvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: FF1, FFF

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL daisie_fill1d(P, 5 * (N + 4 + 2 * kk), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics
    
!  xx1 = c(0,0,x[1:lx],0)
!  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
!  xx3 = x[2 * lx + 1]
!  nil2lx = 3:(lx + 2)
!  il1 = nil2lx+kk-1
!  il2 = nil2lx+kk+1
!  il3 = nil2lx+kk
!  il4 = nil2lx+kk-2
!  in1 = nil2lx+2*kk-1
!  in2 = nil2lx+1
!  in3 = nil2lx+kk
!  ix1 = nil2lx-1
!  ix2 = nil2lx+1
!  ix3 = nil2lx
!  ix4 = nil2lx-2
! new: il3in3 = il3 = in3
! new: in2ix2 = in2 = ix2

      xx1(1) = 0
      xx1(2) = 0
      xx2(1) = 0
	    xx2(2) = 0
      DO I = 3, N + 2
	      xx1(I) = Conc(I - 2)
	      xx2(I) = Conc(N + I - 2)
		    il1(I - 2) = I + kk - 1
		    il2(I - 2) = I + kk + 1
		    il3in3(I - 2) = I + kk
		    il4(I - 2) = I + kk - 2
		    in1(I - 2) = I + 2 * kk - 1
		    in2ix2(I - 2) = I + 1
		    ix1(I - 2) = I - 1
		    ix3(I - 2) = I
		    ix4(I - 2) = I - 2
	    ENDDO
	    xx1(N + 3) = 0
	    xx2(N + 3) = 0
      xx3 = Conc(2 * N + 1)

      DO I = 1, N + 4 + 2 * kk
       laavec(I) = P(I)
	     lacvec(I) = P(I + N + 4 + 2 * kk)
       muvec(I)  = P(I + 2 * (N + 4 + 2 * kk))
	     gamvec(I) = P(I + 3 * (N + 4 + 2 * kk))
       nn(I)     = P(I + 4 * (N + 4 + 2 * kk))
      ENDDO

!  dx1 = laavec[il1 + 1] * xx2[ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3] +
!    lacvec[il1] * nn[in1] * xx1[ix1] + muvec[il2] * nn[in2] * xx1[ix2] +
!    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
!    -gamvec[il3] * xx1[ix3]
!  dx1[1] = dx1[1] + laavec[il3[1]] * xx3 * (kk == 1)
!  dx1[2] = dx1[2] + 2 * lacvec[il3[1]] * xx3 * (kk == 1)
! dx2 = gamvec[il3] * xx1[ix3] +
!    lacvec[il1 + 1] * nn[in1] * xx2[ix1] + muvec[il2 + 1] * nn[in2] * xx2[ix2] +
!    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
!    -laavec[il3 + 1] * xx2[ix3]

      DO I = 1, N 
  	    FF1 = laavec(il1(I) + 1) * xx2(ix1(I))
  	    FF1 = FF1 + lacvec(il4(I) + 1) * xx2(ix4(I))
	      FF1 = FF1 + muvec(il2(I) + 1) * xx2(ix3(I))
	      FF1 = FF1 + lacvec(il1(I)) * nn(in1(I)) * xx1(ix1(I))
	      FF1 = FF1 + muvec(il2(I)) * nn(in2ix2(I)) * xx1(in2ix2(I))
	      FFF = (muvec(il3in3(I)) + lacvec(il3in3(I)))
	      FF1 = FF1 - FFF * nn(il3in3(I)) * xx1(ix3(I))
	      FF1 = FF1 - gamvec(il3in3(I)) * xx1(ix3(I))
	      dConc(I) = FF1
        FF1 = gamvec(il3in3(I)) * xx1(ix3(I))
		    FF1 = FF1 + lacvec(il1(I) + 1) * nn(in1(I)) * xx2(ix1(I))
		    FF1 = FF1 + muvec(il2(I)+1) * nn(in2ix2(I)) * xx2(in2ix2(I))
		    FFF = muvec(il3in3(I) + 1) + lacvec(il3in3(I) + 1)
		    FF1 = FF1 - FFF * nn(il3in3(I) + 1) * xx2(ix3(I))
		    FF1 = FF1 - laavec(il3in3(I) + 1) * xx2(ix3(I))
        dConc(N + I) = FF1
      ENDDO
	    IF(kk .EQ. 1) THEN
	      dConc(1) = dConc(1) + laavec(il3in3(1)) * xx3
	      dConc(2) = dConc(2) + 2 * lacvec(il3in3(1)) * xx3
      ENDIF

!  dx3 = -(laavec[il3in3[1]] + lacvec[il3in3[1]] + gamvec[il3in3[1]] + muvec[il3in3[1]]) * xx3

      FFF = laavec(il3in3(1)) + lacvec(il3in3(1))
	    FFF = FFF + gamvec(il3in3(1)) + muvec(il3in3(1))
	    dConc(2 * N + 1) = -1 * FFF * xx3

      END SUBROUTINE daisie_runmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE daisie_runmod2 (neq, t, Conc, dConc, yout, ip)
      USE daisie_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(3 * N), dConc(3 * N), yout(*)
      DOUBLE PRECISION  :: xx1(N + 3), xx2(N + 3), xx3(N + 3)
	    INTEGER 		      :: il1(N), il2(N), il3in3(N), il4(N)
	    INTEGER           :: in1(N), in2ix2(N), in4ix1(N)
      INTEGER           :: ix3(N), ix4(N)
      DOUBLE PRECISION  :: laavec(N + 4 + 2 * kk),lacvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: muvec(N + 4 + 2 * kk),gamvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: FF1, FFF

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL daisie_fill1d(P, 5 * (N + 4 + 2 * kk), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics
    
!  xx1 = c(0,0,x[1:lx],0)
!  xx2 = c(0,0,x[(lx + 1):(2 * lx)],0)
!  xx3 = c(0,0,x[(2 * lx + 1):(3 * lx)],0)
!  nil2lx = 3:(lx + 2)
!  il1 = nil2lx+kk-1
!  il2 = nil2lx+kk+1
!  il3 = nil2lx+kk
!  il4 = nil2lx+kk-2
!  in1 = nil2lx+2*kk-1
!  in2 = nil2lx+1
!  in3 = nil2lx+kk
!  in4 = nil2lx-1
!  ix1 = nil2lx-1
!  ix2 = nil2lx+1
!  ix3 = nil2lx
!  ix4 = nil2lx-2
! new: il3in3 = il3 = in3
! new: in2ix2 = in2 = ix2
! new: in4ix1 = in4 = ix1

      xx1(1) = 0
      xx1(2) = 0
      xx2(1) = 0
	    xx2(2) = 0
      xx3(1) = 0
	    xx3(2) = 0
      DO I = 3, N + 2
	      xx1(I) = Conc(I - 2)
	      xx2(I) = Conc(N + I - 2)
	      xx3(I) = Conc(2 * N + I - 2)
		    il1(I - 2) = I + kk - 1
		    il2(I - 2) = I + kk + 1
		    il3in3(I - 2) = I + kk
		    il4(I - 2) = I + kk - 2
		    in1(I - 2) = I + 2 * kk - 1
		    in2ix2(I - 2) = I + 1
		    in4ix1(I - 2) = I - 1
		    ix3(I - 2) = I
		    ix4(I - 2) = I - 2
	    ENDDO
	    xx1(N + 3) = 0
	    xx2(N + 3) = 0
      xx3(N + 3) = 0

      DO I = 1, N + 4 + 2 * kk
       laavec(I) = P(I)
	     lacvec(I) = P(I + N + 4 + 2 * kk)
       muvec(I)  = P(I + 2 * (N + 4 + 2 * kk))
	     gamvec(I) = P(I + 3 * (N + 4 + 2 * kk))
       nn(I)     = P(I + 4 * (N + 4 + 2 * kk))
      ENDDO

! dx1 = (laavec[il3] * xx3[ix3] + 2 * lacvec[il1] * xx3[in4ix1]) * (kk == 1) + 
!    laavec[il1 + 1] * xx2[in4ix1] + lacvec[il4 + 1] * xx2[ix4] + muvec[il2 + 1] * xx2[ix3]
!    lacvec[il1] * nn[in1] * xx1[in4ix1] + muvec[il2] * nn[in2] * xx1[in2ix2] +
!    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] - gamvec[il3] * xx1[ix3]
! dx2 = gamvec[il3] * xx1[ix3] +
!    lacvec[il1 + 1] * nn[in1] * xx2[in4ix1] + muvec[il2 + 1] * nn[in2] * xx2[in2ix2] +
!    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
!    -laavec[il3 + 1] * xx2[ix3]
! dx3 = lacvec[il1] * nn[in4ix1] * xx3[in4ix1] + muvec[il2] * nn[in2] * xx3[in2ix2] +
!    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
!    -(laavec[il3] + gamvec[il3]) * xx3[ix3]
  
      DO I = 1, N 
  	    FF1 = laavec(il1(I) + 1) * xx2(in4ix1(I))
  	    FF1 = FF1 + lacvec(il4(I) + 1) * xx2(ix4(I))
	      FF1 = FF1 + muvec(il2(I) + 1) * xx2(ix3(I))
	      FF1 = FF1 + lacvec(il1(I)) * nn(in1(I)) * xx1(in4ix1(I))
	      FF1 = FF1 + muvec(il2(I)) * nn(in2ix2(I)) * xx1(in2ix2(I))
	      FFF = muvec(il3in3(I)) + lacvec(il3in3(I))
	      FF1 = FF1 - FFF * nn(il3in3(I)) * xx1(ix3(I))
	      FF1 = FF1 - gamvec(il3in3(I)) * xx1(ix3(I))
	      FFF = 0	      
	      IF(kk .EQ. 1) THEN
	         FFF = laavec(il3in3(I)) * xx3(ix3(I))
 	         FFF = FFF + 2 * lacvec(il1(I)) * xx3(in4ix1(I))
	      ENDIF
	      dConc(I) = FF1 + FFF
        FF1 = gamvec(il3in3(I)) * xx1(ix3(I))
		    FF1 = FF1 + lacvec(il1(I) + 1) * nn(in1(I)) * xx2(in4ix1(I))
		    FF1 = FF1 + muvec(il2(I)+1) * nn(in2ix2(I)) * xx2(in2ix2(I))
		    FFF = muvec(il3in3(I) + 1) + lacvec(il3in3(I) + 1)
		    FF1 = FF1 - FFF * nn(il3in3(I) + 1) * xx2(ix3(I))
		    FF1 = FF1 - laavec(il3in3(I) + 1) * xx2(ix3(I))
        dConc(N + I) = FF1
        FF1 = lacvec(il1(I)) * nn(in4ix1(I)) * xx3(in4ix1(I))
        FF1 = FF1 + muvec(il2(I)) * nn(in2ix2(I)) * xx3(in2ix2(I))
        FFF = lacvec(il3in3(I)) + muvec(il3in3(I))
        FF1 = FF1 - FFF * nn(il3in3(I)) * xx3(ix3(I))
        FF1 = FF1-(laavec(il3in3(I))+gamvec(il3in3(I)))*xx3(ix3(I))
        dConc(2 * N + I) = FF1
      ENDDO

      END SUBROUTINE daisie_runmod2
      
