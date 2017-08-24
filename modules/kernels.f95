MODULE kernels

! Module to calculate geometric sources
USE geometric, ONLY: polyprism_class


IMPLICIT NONE
! Variables:
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(p=4, r=4)    ! Single precision 
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(p=8, r=8)    ! Double precision                   
INTEGER(KIND=DP):: i, j, k


CONTAINS

! CREATE THE PLOUFF CALCULATIONS V1, V2, V3, V4, V5, V6

! SUBROUTINE KERNELL1(V1)
! SUBROUTINE kERNEL2 (V2) ......


END MODULE kernels
