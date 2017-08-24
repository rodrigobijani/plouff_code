MODULE geometric
! Module to generate geometric sources, prisms
IMPLICIT NONE
! Variables:
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(p=4, r=4)    ! Single precision 
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(p=8, r=8)    ! Double precision                   
INTEGER(KIND=DP):: i, j, k

CONTAINS

! Create a class with geometric info about the polygon:
TYPE polyprism_class
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE:: xv,yv   ! vertices of the prism (horizontal coordinates x,y)   
 REAL(KIND=DP):: z1,z2 ! top and bottom of the prism
 INTEGER:: nv ! number of vertices
 CHARACTER(LEN=4):: props
END TYPE polyprism_class

! *********** CLEAN AND DEALLOCATE OBJECT *************************! 
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE geometric_clear(pp)
 ! Clears the object.
 IMPLICIT NONE
 TYPE(polyprism_class), INTENT(INOUT) :: pp

 ! Clear all allocated arrays:
 CALL geometric_deallocate(pp)
 ! Clear everything else:
 pp%z1    = 0.0
 pp%z2    = 0.0
 pp%nv    = 0
 pp%props = 'none'

END SUBROUTINE geometric_clear

SUBROUTINE geometric_deallocate(pp)
! Deallocates the arrays in the object. 
 IMPLICIT NONE
 TYPE(polyprism_class), INTENT(INOUT):: pp
 IF (ALLOCATED(pp%xv))  DEALLOCATE(pp%xv)
 IF (ALLOCATED(pp%yv))  DEALLOCATE(pp%yv)  
END SUBROUTINE geometric_deallocate


END MODULE geometric
