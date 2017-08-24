MODULE gravmag
! Module to generate geometric sources
IMPLICIT NONE
! Variables:
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(p=8, r=8)  ! DP = Double Precision
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(p=4, r=4)  ! SP = Single Precision
INTEGER(KIND=DBL):: i, j, k


CONTAINS


SUBROUTINE mag_plouff(...)
! computes the magnetic field of a polygonal prism in three dimensions:
IMPLICIT NONE


! Aqui calculamos os 6 Vs do problema:

!V1 = DO i=1,N (S(x1,x2,y1,y2)CF - C²W) ENDDO
             



CONTAINS

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
FUNCTION S(x1,x2,y1,y2)
IMPLICIT NONE
REAL(KIND=DP):: x1,x2,y1,y2,S

S = (x2-x1) / DSQRT( (x2-x1)**2 + (y2-y1)**2 )
! Check denominator condition:
IF( (x2-x1) .AND. (y2-y1) == 0.0)THEN
  PRINT*,'FUNCTION S: - COORDINATE ERROR: CHECK INPUT FILE'
  STOP
ENDIF
END FUNCTION S

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
FUNCTION C(x1,x2,y1,y2)
IMPLICIT NONE
REAL(KIND=DP):: x1,x2,y1,y2,C

C = (y2-y1) / DSQRT( (x2-x1)**2 + (y2-y1)**2 )
! Check denominator condition:
IF( (x2-x1) .AND. (y2-y1) == 0.0)THEN
  PRINT*,'FUNCTION C: - COORDINATE ERROR: CHECK INPUT FILE'
  STOP
ENDIF
END FUNCTION C

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
FUNCTION F(z1,z2,R11,R12,R21,R22)
IMPLICIT NONE
REAL(KIND=DP):: R11,R12,R21,R22,z1,z2,F

F = DLOG( (R22 + z2 * R11 + z1) / (R12 + z2 * R21 + z1) )
END FUNCTION F

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
FUNCTION P(x1,x2,y1,y2)
IMPLICIT NONE
REAL(KIND=DP):: x1,x2,y1,y2

P = (x1*y2 - x2*y1) / DSQRT( (x2-y1)**2 + (y2-y1)**2 )
! Check denominator condition:
IF( (x2-x1) .AND. (y2-y1) == 0.0)THEN
  PRINT*,'FUNCTION P: - COORDINATE ERROR: CHECK INPUT FILE'
  STOP
ENDIF

END FUNCTION P
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

! FALTA W 



END SUBROUTINE mag_plouff




END MODULE gravmag
