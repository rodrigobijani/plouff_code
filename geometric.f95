MODULE geometric
! Module to generate geometric sources
IMPLICIT NONE
! Variables:
INTEGER, PARAMETER:: SGL = SELECTED_REAL_KIND(p=8, r=8)
INTEGER, PARAMETER:: DBL = SELECTED_REAL_KIND(p=4, r=4)
INTEGER(KIND=DBL):: i, j, k


CONTAINS

! Create a class with geometric info about the polygon:
!TYPE coords
!   PUBLIC
!   REAL(KIND=DP), DIMENSION(2):: x,y,z ! coordinates of a 3D polygon  
!END TYPE coords




END MODULE geometric
