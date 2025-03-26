!
!  Copyright (c) 2023, Rampal S. Etienne
!
!  Distributed under the Boost Software License, Version 1.0. (See
!  accompanying file LICENSE_1_0.txt or copy at
!  http://www.boost.org/LICENSE_1_0.txt)
!
!==========================================================================
! Helper function:
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE daisie_fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
      II = II
      vec(1:DIMP) = parms(II + 1:II + DIMP)

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
      !COMMON /XCBPar/parms                 ! common block

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
      INTEGER           :: il1(N), il2(N), il3in3(N), il4(N)
      INTEGER           :: in1(N), in2ix2(N)
      INTEGER           :: ix1(N), ix3(N), ix4(N)
      DOUBLE PRECISION  :: laavec(N + 4 + 2 * kk),lacvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: muvec(N + 4 + 2 * kk),gamvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 4 + 2 * kk)

! parameters - named here
      !DOUBLE PRECISION rn(2)
      !COMMON /XCBPar/rn

! local variables
      !CHARACTER(len=100) msg

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

!  dx1 = laavec[il1 + 1] * xx2[ix1] +
!   lacvec[il4 + 1] * xx2[ix4] +
!    muvec[il2 + 1] * xx2[ix3] +
!    lacvec[il1] * nn[in1] * xx1[ix1] +
!    muvec[il2] * nn[in2] * xx1[ix2] +
!    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
!    -gamvec[il3] * xx1[ix3]
! dx2 = gamvec[il3] * xx1[ix3] +
!    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
!    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
!    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
!    -laavec[il3 + 1] * xx2[ix3]
! dx3 <- 0

      DO I = 1, N
            dConc(I) = &
                laavec(il1(I) + 1) * xx2(ix1(I)) + &
                lacvec(il4(I) + 1) * xx2(ix4(I)) + &
                muvec(il2(I) + 1) * xx2(ix3(I)) + &
                lacvec(il1(I)) * nn(in1(I)) * xx1(ix1(I)) + &
                muvec(il2(I)) * nn(in2ix2(I)) * xx1(in2ix2(I)) - &
                (muvec(il3in3(I)) + lacvec(il3in3(I))) * &
                nn(il3in3(I)) * xx1(ix3(I)) - &
                gamvec(il3in3(I)) * xx1(ix3(I))
        dConc(N + I) = &
            gamvec(il3in3(I)) * xx1(ix3(I)) + &
            lacvec(il1(I) + 1) * nn(in1(I)) * xx2(ix1(I)) + &
            muvec(il2(I) + 1) * nn(in2ix2(I)) * xx2(in2ix2(I)) - &
            (muvec(il3in3(I) + 1) + lacvec(il3in3(I) + 1)) * &
            nn(il3in3(I) + 1) * xx2(ix3(I)) - &
            laavec(il3in3(I) + 1) * xx2(ix3(I))
            dConc(2 * N + 1) = 0

      ENDDO

      END SUBROUTINE daisie_runmod

!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================

      SUBROUTINE daisie_runmod1 (neq, t, Conc, dConc, yout, ip)
      USE daisie_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(4 * N), dConc(4 * N), yout(*)
      DOUBLE PRECISION  :: xx1(N + 3), xx2(N + 3), xx3(N + 3), xx4(N + 3)
      INTEGER           :: il1(N), il2(N), il3in3(N), il4(N)
      INTEGER           :: in1(N), in2ix2(N), in4ix1(N)
      INTEGER           :: ix3(N), ix4(N)
      DOUBLE PRECISION  :: laavec(N + 4 + 2 * kk),lacvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: muvec(N + 4 + 2 * kk),gamvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: FF1, FFF

! parameters - named here
      DOUBLE PRECISION rn(2)
      !COMMON /XCBPar/rn

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

!  xx1 <- c(0,0,x[1:lx],0)
!  xx2 <- c(0,0,x[(lx + 1):(2 * lx)],0)
!  xx3 <- c(0,0,x[(2 * lx + 1):(3 * lx)],0)
!  xx4 <- c(0,0,x[(3 * lx + 1):(4 * lx)],0)
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
      xx4(1) = 0
      xx4(2) = 0
      DO I = 3, N + 2
            xx1(I) = Conc(I - 2)
            xx2(I) = Conc(N + I - 2)
            xx3(I) = Conc(2 * N + I - 2)
            xx4(I) = Conc(3 * N + I - 2)
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
      xx4(N + 3) = 0

      DO I = 1, N + 4 + 2 * kk
            laavec(I) = P(I)
            lacvec(I) = P(I + N + 4 + 2 * kk)
            muvec(I)  = P(I + 2 * (N + 4 + 2 * kk))
            gamvec(I) = P(I + 3 * (N + 4 + 2 * kk))
            nn(I)     = P(I + 4 * (N + 4 + 2 * kk))
      ENDDO

!  dx1 <- lacvec[il1] * xx1[ix1] +
!    laavec[il1 + 1] * xx2[ix1] +
!    lacvec[il4 + 1] * xx2[ix4] +
!    muvec[il2] * nn[in2] * xx1[ix2] +
!    muvec[il3 + 1] * xx2[ix3] +
!    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
!    -gamvec[il3] * xx1[ix3]

      DO I = 1, N
            dConc(I) = &
               lacvec(il1(I)) * nn(in1(I)) * xx1(in4ix1(I)) + &
               laavec(il1(I) + 1) * xx2(in4ix1(I)) + &
               lacvec(il4(I) + 1) * xx2(ix4(I)) + &
               muvec(il2(I)) * nn(in2ix2(I)) * xx1(in2ix2(I)) + &
               muvec(il3in3(I) + 1) * xx2(ix3(I)) - &
               (muvec(il3in3(I)) + lacvec(il3in3(I))) * &
               nn(il3in3(I)) * xx1(ix3(I)) - &
               gamvec(il3in3(I)) * xx1(ix3(I))

!  dx2 <- gamvec[il3] * xx1[ix3] +
!    gamvec[il3] * xx3[ix3] +
!    gamvec[il3 + 1] * xx4[ix3] +
!    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
!    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
!    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
!    -laavec[il3 + 1] * xx2[ix3]

        dConc(N + I) = &
           gamvec(il3in3(I)) * xx1(ix3(I)) + &
           gamvec(il3in3(I)) * xx3(ix3(I)) + &
           gamvec(il3in3(I) + 1) * xx4(ix3(I)) + &
                   lacvec(il1(I) + 1) * nn(in1(I)) * xx2(in4ix1(I)) + &
                   muvec(il2(I) + 1) * nn(in2ix2(I)) * xx2(in2ix2(I)) - &
                   (muvec(il3in3(I) + 1) + lacvec(il3in3(I) + 1)) * &
                   nn(il3in3(I) + 1) * xx2(ix3(I)) - &
                   laavec(il3in3(I) + 1) * xx2(ix3(I))

!  dx3 <- lacvec[il1] * nn[in1] * xx3[ix1] +
!    laavec[il1 + 1] * xx4[ix1] +
!    lacvec[il4 + 1] * xx4[ix4] +
!    muvec[il2] * nn[in2] * xx3[ix2] +
!    muvec[il3 + 1] * xx4[ix3] +
!    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
!    -gamvec[il3] * xx3[ix3]

        dConc(2 * N + I) = &
          lacvec(il1(I)) * nn(in1(I)) * xx3(in4ix1(I)) + &
          laavec(il1(I) + 1) * xx4(in4ix1(I)) + &
          lacvec(il4(I) + 1) * xx4(ix4(I)) + &
          muvec(il2(I)) * nn(in2ix2(I)) * xx3(in2ix2(I)) + &
          muvec(il3in3(I) + 1) * xx4(ix3(I)) - &
          (lacvec(il3in3(I)) + muvec(il3in3(I))) * &
          nn(il3in3(I)) * xx3(ix3(I)) - &
          gamvec(il3in3(I)) * xx3(ix3(I))

!  dx4 <- lacvec[il1 + 1] * nn[in1] * xx4[ix1] +
!    muvec[il2 + 1] * nn[in2] * xx4[ix2] +
!    -(laavec[il3 + 1]) + lacvec[il3 + 1] + muvec[il3 + 1]) * &
!    nn[in3 + 1] * xx4[ix3] +
!    -gamvec[il3 + 1] * xx4[ix3]

      dConc(3 * N + I) = &
           lacvec(il1(I) + 1) * nn(in1(I)) * xx4(in4ix1(I)) + &
           muvec(il2(I) + 1) * nn(in2ix2(I)) * xx4(in2ix2(I)) - &
           (lacvec(il3in3(I) + 1) + muvec(il3in3(I) + 1)) * &
           nn(il3in3(I) + 1) * xx4(ix3(I)) - &
           laavec(il3in3(I) + 1) * xx4(ix3(I)) - &
           gamvec(il3in3(I) + 1) * xx4(ix3(I))

      ENDDO

    END SUBROUTINE daisie_runmod1


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
      INTEGER           :: il1(N), il2(N), il3in3(N), il4(N)
      INTEGER           :: in1(N), in2ix2(N), in4ix1(N)
      INTEGER           :: ix3(N), ix4(N)
      DOUBLE PRECISION  :: laavec(N + 4 + 2 * kk),lacvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: muvec(N + 4 + 2 * kk),gamvec(N + 4 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 4 + 2 * kk)

! parameters - named here
      !DOUBLE PRECISION rn(2)
      !COMMON /XCBPar/rn

! local variables
      !CHARACTER(len=100) msg

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

!  dx1 = (laavec[il3] * xx3[ix3] +
!    2 * lacvec[il1] * xx3[ix1]) * (kk == 1) +
!    laavec[il1 + 1] * xx2[ix1] +
!    lacvec[il4 + 1] * xx2[ix4] +
!    muvec[il2 + 1] * xx2[ix3] +
!    lacvec[il1] * nn[in1] * xx1[ix1] +
!    muvec[il2] * nn[in2] * xx1[ix2] +
!    -(muvec[il3] + lacvec[il3]) * nn[in3] * xx1[ix3] +
!    -gamvec[il3] * xx1[ix3]

!  dx2 <- gamvec[il3] * xx1[ix3] +
!    lacvec[il1 + 1] * nn[in1] * xx2[ix1] +
!    muvec[il2 + 1] * nn[in2] * xx2[ix2] +
!    -(muvec[il3 + 1] + lacvec[il3 + 1]) * nn[in3 + 1] * xx2[ix3] +
!    -laavec[il3 + 1] * xx2[ix3]

!  dx3 <- lacvec[il1] * nn[in4] * xx3[ix1] +
!    muvec[il2] * nn[in2] * xx3[ix2] +
!    -(lacvec[il3] + muvec[il3]) * nn[in3] * xx3[ix3] +
!    -(laavec[il3] + gamvec[il3]) * xx3[ix3]

      DO I = 1, N
            dConc(I) = 0
            IF(kk .EQ. 1) THEN
               dConc(I) = laavec(il3in3(I)) * xx3(ix3(I)) + &
                2 * lacvec(il1(I)) * xx3(in4ix1(I))
            ENDIF
            dConc(I) = dConc(I) + &
                laavec(il1(I) + 1) * xx2(in4ix1(I)) + &
                lacvec(il4(I) + 1) * xx2(ix4(I)) + &
                muvec(il2(I) + 1) * xx2(ix3(I)) + &
                lacvec(il1(I)) * nn(in1(I)) * xx1(in4ix1(I)) + &
                muvec(il2(I)) * nn(in2ix2(I)) * xx1(in2ix2(I)) - &
                (muvec(il3in3(I)) + lacvec(il3in3(I))) * &
                nn(il3in3(I)) * xx1(ix3(I)) - &
                gamvec(il3in3(I)) * xx1(ix3(I))
            dConc(N + I) = &
                  gamvec(il3in3(I)) * xx1(ix3(I)) + &
                  lacvec(il1(I) + 1) * nn(in1(I)) * xx2(in4ix1(I)) + &
                  muvec(il2(I) + 1) * nn(in2ix2(I)) * xx2(in2ix2(I)) - &
                  (muvec(il3in3(I) + 1) + lacvec(il3in3(I) + 1)) * &
                  nn(il3in3(I) + 1) * xx2(ix3(I)) - &
                  laavec(il3in3(I) + 1) * xx2(ix3(I))
            dConc(2 * N + I) = &
                  lacvec(il1(I)) * nn(in4ix1(I)) * xx3(in4ix1(I)) + &
                  muvec(il2(I)) * nn(in2ix2(I)) * xx3(in2ix2(I)) - &
                  (lacvec(il3in3(I)) + muvec(il3in3(I))) * &
                  nn(il3in3(I)) * xx3(ix3(I)) - &
                  (laavec(il3in3(I)) + gamvec(il3in3(I))) * xx3(ix3(I))
      ENDDO

      END SUBROUTINE daisie_runmod2

!==========================================================================
! module with declarations
!==========================================================================

      MODULE daisie_dimmod3

      ! length of the vector -  decided in R-code
      INTEGER  :: N1
      INTEGER  :: N2

      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)

      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE daisie_dimmod3


!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE daisie_initmod3 (steadyparms)
      USE daisie_dimmod3

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 2  ! constant-length parameters

      DOUBLE PRECISION parms(nparsmall)
      !COMMON /XCBPar/parms                 ! common block

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector
      N1 = INT(parms(1) + 1e-6)
      N2 = INT(parms(2) + 1e-6)

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)
      ALLOCATE(P(7 * N1 + 14 * N1 * N2 + 2))
! 2 for N1 and N2
! 3 * N for a1-a3
! 8 * N^2 + 3 * N for b1-b11
! 6 * N^2 + N for c1-c7
! Total is 7 * N + 14 * N^2 + 2
      initialised = .FALSE.

      END SUBROUTINE daisie_initmod3


!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================

      SUBROUTINE daisie_runmod3 (neq, t, Conc, dConc, yout, ip)
      USE daisie_dimmod3
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii, j
      DOUBLE PRECISION  :: t, yout(*)
      DOUBLE PRECISION  :: Conc(N1 + 2 * N1 * N2), dConc(N1 + 2 * N1 * N2)
      DOUBLE PRECISION  :: xx1(N1 + 3), xx2(N1 + 3,N2 + 3), xx3(N1 + 3, N2 + 3)
      INTEGER           :: nil2lx1(N1), nil2lx2(N2)
      DOUBLE PRECISION  :: a1(N1), a2(N1), a3(N1), b4(N1), b5(N1), b6(N1),c2(N1)
      DOUBLE PRECISION  :: b1(N1, N2), b2(N1, N2), b3(N1, N2), b7(N1, N2)
      DOUBLE PRECISION  :: b8(N1, N2), b9(N1, N2), b10(N1, N2), b11(N1, N2)
      DOUBLE PRECISION  :: c1(N1, N2), c3(N1, N2), c4(N1, N2), c5(N1, N2)
      DOUBLE PRECISION  :: c6(N1, N2), c7(N1, N2)
! parameters - named here
      !DOUBLE PRECISION rn(2)
      !COMMON /XCBPar/rn

! local variables
      !CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough")

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL daisie_fill1d(P, 7 * N1 + 14 * N1 * N2, yout, ii) ! ii is updated in fill1d
        Initialised = .TRUE.      ! to prevent from initialising more than once
      ENDIF

! dynamics

!  lx1 <- cp$lx1
!  lx2 <- cp$lx2
!  x1 <- x[1:lx1]
!  x2 <- x[(lx1 + 1):(lx1 + lx1 * lx2)]
!  x3 <- x[(lx1 + lx1 * lx2 + 1):(lx1 + 2 * lx1 * lx2)]
!  dim(x2) <- c(lx1 ,lx2)
!  dim(x3) <- c(lx1, lx2)
!  xx1 <- rep(0,lx1 + 3)
!  xx2 <- array(0,dim = c(lx1 + 3, lx2 + 3))
!  xx3 <- array(0,dim = c(lx1 + 3, lx2 + 3))
!  nil2lx1 <- 3:(lx1 + 2)
!  nil2lx2 <- 3:(lx2 + 2)
!  xx1[nil2lx1] <- x1
!  xx2[nil2lx1, nil2lx2] <- x2
!  xx3[nil2lx1, nil2lx2] <- x3
!  dx1 <-
!    cp$a1 * xx1[nil2lx1 - 1] +
!    cp$a2 * xx1[nil2lx1 + 1] -
!    cp$a3 * xx1[nil2lx1]
!  dx2 <-
!    cp$b1 * xx3[nil2lx1, nil2lx2 - 1] +
!    cp$b2 * xx3[nil2lx1, nil2lx2 - 2] +
!    cp$b3 * xx3[nil2lx1, nil2lx2] +
!    cp$b7 * xx2[nil2lx1 - 1, nil2lx2]  +
!    cp$b8 * xx2[nil2lx1, nil2lx2 - 1] +
!    cp$b9 * xx2[nil2lx1 + 1, nil2lx2] +
!    cp$b10 * xx2[nil2lx1, nil2lx2 + 1] -
!    cp$b11 * xx2[nil2lx1, nil2lx2]
!  dx2[,1] <- dx2[,1] +
!    cp$b4 * xx1[nil2lx1] +
!    cp$b5 * xx1[nil2lx1 - 1] +
!    cp$b6 * xx1[nil2lx1 - 2]
!  dim(dx2) <- c(lx1 * lx2, 1)
!  dx3 <-
!    cp$c1 * xx2[nil2lx1, nil2lx2] +
!    cp$c3 * xx3[nil2lx1 - 1, nil2lx2] +
!    cp$c4 * xx3[nil2lx1, nil2lx2 - 1]  +
!    cp$c5 * xx3[nil2lx1 + 1, nil2lx2] +
!    cp$c6 * xx3[nil2lx1, nil2lx2 + 1] -
!    cp$c7 * xx3[nil2lx1, nil2lx2]
!  dx3[,1] <- dx3[,1] +
!    cp$c2 * xx1[nil2lx1]
!  dim(dx3) <- c(lx1 * lx2, 1)

      xx1 = 0.0
      xx2 = 0.0
      xx3 = 0.0
      nil2lx1 = 3:(N1 + 2)
      nil2lx2 = 3:(N2 + 2)
      xx1(nil2lx1 = Conc(1:N1)
      xx2(nil2lx1), nil2lx2) = &
         RESHAPE(Conc(N1 + 1:N1 + N1 * N2), [N1, N2])
      xx3(nil2lx1), nil2lx2) = &
         RESHAPE(Conc(N1 + N1 * N2 + 1:N1 + 2 * N1 * N2), [N1, N2])

      a1 = P(1:N1)
      a2 = P(N1 + 1:2 * N1)
      a3 = P(2 * N1 + 1:3 * N1)
      b1 = RESHAPE(P(3 * N1 + 1:3 * N1 + N1 * N2),[N1, N2])
      b2 = RESHAPE(P(3 * N1 + N1 * N2 + 1:3 * N1 + 2 * N1 * N2),[N1, N2])
      b3 = RESHAPE(P(3 * N1 + 2 * N1 * N2 + 1:3 * N1 + 3 * N1 * N2),[N1, N2])
      b4 = P(3 * N1 + 3 * N1 * N2 + 1:4 * N1 + 3 * N1 * N2)
      b5 = P(4 * N1 + 3 * N1 * N2 + 1:5 * N1 + 3 * N1 * N2)
      b6 = P(5 * N1 + 3 * N1 * N2 + 1:6 * N1 + 3 * N1 * N2)
      b7 = RESHAPE(P(6 * N1 + 3 * N1 * N2 + 1:6 * N1 + 4 * N1 * N2), [N1, N2])
      b8 = RESHAPE(P(6 * N1 + 4 * N1 * N2 + 1:6 * N1 + 5 * N1 * N2), [N1, N2])
      b9 = RESHAPE(P(6 * N1 + 5 * N1 * N2 + 1:6 * N1 + 6 * N1 * N2), [N1, N2])
      b10 = RESHAPE(P(6 * N1 + 6 * N1 * N2 + 1:6 * N1 + 7 * N1 * N2), [N1, N2])
      b11 = RESHAPE(P(6 * N1 + 7 * N1 * N2 + 1:6 * N1 + 8 * N1 * N2), [N1, N2])
      c1 = RESHAPE(P(6 * N1 + 8 * N1 * N2 + 1:6 * N1 + 9 * N1 * N2), [N1, N2])
      c2 = P(6 * N1 + 9 * N1 * N2 + 1:7 * N1 + 9 * N1 * N2)
      c3 = RESHAPE(P(7 * N1 + 9 * N1 * N2 + 1:7 * N1 + 10 * N1 * N2), [N1, N2])
      c4 = RESHAPE(P(7 * N1 + 10 * N1 * N2 + 1:7 * N1 + 11 * N1 * N2), [N1, N2])
      c5 = RESHAPE(P(7 * N1 + 11 * N1 * N2 + 1:7 * N1 + 12 * N1 * N2), [N1, N2])
      c6 = RESHAPE(P(7 * N1 + 12 * N1 * N2 + 1:7 * N1 + 13 * N1 * N2), [N1, N2])
      c7 = RESHAPE(P(7 * N1 + 13 * N1 * N2 + 1:7 * N1 + 14 * N1 * N2), [N1, N2])

      dConc(1:N1) = a1 * xx1(nil2lx1 - 1) + &
                    a2 * xx1(nil2lx1 + 1) - &
                    a3 * xx1(nil2lx1)
      dConc(N1 + 1:N1 + 2 * N1 * N2) = b1 * xx3(nil2lx1, nil2lx2 - 1) + &
                                       b2 * xx3(nil2lx1, nil2lx2 - 2) + &
                                       b3 * xx3(nil2lx1, nil2lx2 + &
                                       b7 * xx2(nil2lx1 - 1, nil2lx2) + &
                                       b8 * xx2(nil2lx1, nil2lx2 - 1) + &
                                       b9 * xx2(nil2lx1 + 1, nil2lx2) + &
                                       b10 * xx2(nil2lx1, nil2lx2 + 1) - &
                                       b11 * xx2(nil2lx1, nil2lx2)
      dConc(N1 + 1:2 * N1) = dConc(N1 + 1:2 * N1) + &
                             b4 * xx1(nil2lx1) + &
                             b5 * xx1(nil2lx1 - 1) + &
                             b6 * xx1(nil2lx1 - 2)
      dConc(N1 + N1 * N2 + 1:N1 + 2 * N1 * N2) = &
                             c1 * xx2(nil2lx1, nil2lx2) + &
                             c3 * xx3(nil2lx1 - 1, nil2lx2) + &
                             c4 * xx3(nil2lx1, nil2lx2 - 1) + &
                             c5 * xx3(nil2lx1 + 1, nil2lx2) + &
                             c6 * xx3(nil2lx1, nil2lx2 + 1) - &
                             c7 * xx3(nil2lx1, nil2lx2)
      dConc(N1 + N1 * N2 + 1:2 * N1 + N1 * N2) = &
                             dConc(N1 + N1 * N2 + 1:2 * N1 + N1 * N2) + &
                             c2 * xx1(nil2lx1)

      END SUBROUTINE daisie_runmod3
