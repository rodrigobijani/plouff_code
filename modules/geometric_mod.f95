MODULE geometric_mod
! Module to generate geometric elements for computing grav and mag data
IMPLICIT NONE
! Variables:
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(p=4, r=4)    ! Single precision 
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(p=8, r=8)    ! Double precision                   
INTEGER(KIND=DP):: i, j, k

! Create a class with geometric info about the polygon:
TYPE polyprism_type
 REAL(KIND=DP), DIMENSION(:), ALLOCATABLE:: xv,yv   ! vertices of the prism (horizontal coordinates x,y)   
 REAL(KIND=DP):: z1,z2 ! top and bottom of the prism
 INTEGER:: nv ! number of vertices
 CHARACTER(LEN=4):: props ! physical property of the prism (sus = susceptibility , dens = density)
 REAL(KIND=DP):: ppval, inc, dec ! value of the physical property  (inc, dec =  magnetization of the polyprism)
END TYPE polyprism_type


CONTAINS

! *********** CLEAN AND DEALLOCATE OBJECT *************************! 
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE geometric_clear(pp)
 ! Clears the object.
 IMPLICIT NONE
 TYPE(polyprism_type), INTENT(INOUT) :: pp

 ! Clear all allocated arrays:
 CALL geometric_clear(pp)
 ! Set default values for all variables in pp object:
 pp%z1    = 0.0
 pp%z2    = 0.0
 pp%nv    = 0
 pp%ppval = 0.0
 pp%props = 'none'
 pp%inc    = 0.0
 pp%dec    = 0.0
END SUBROUTINE geometric_clear
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE geometric_deallocate(pp)
! Deallocates the arrays in the object. 
 IMPLICIT NONE
 TYPE(polyprism_type), INTENT(INOUT):: pp
 IF (ALLOCATED(pp%xv))  DEALLOCATE(pp%xv)
 IF (ALLOCATED(pp%yv))  DEALLOCATE(pp%yv)  
END SUBROUTINE geometric_deallocate
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE geometric_allocate(pp)
IMPLICIT NONE
TYPE(polyprism_type), INTENT(INOUT):: pp
INTEGER(KIND=DP):: nv,ierr
! Allocate the arrays in pp:
nv = pp%nv
ALLOCATE(pp%xv(nv), pp%yv(nv), STAT=ierr )
IF(ierr /= 0)THEN
  PRINT*, 'error while allocating arrays of polyprism_type'
  RETURN
ENDIF
! initializing arrays:
pp%xv = 0.0
pp%yv = 0.0
END SUBROUTINE geometric_allocate
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

SUBROUTINE geometric_read_inputs(pp, fid, inputfile)
! reads input file with info about the polyprism 
IMPLICIT NONE
TYPE(polyprism_type), INTENT(INOUT):: pp  ! polyprism object
CHARACTER(LEN=*), INTENT(IN):: inputfile ! name of inputfile
INTEGER(KIND=DP), INTENT(IN):: fid ! number of the inputfile
LOGICAL:: ok
INTEGER:: i, ierr

! CHECK IF THE FILE REALLY EXISTIS:
INQUIRE(FILE=TRIM(inputfile),EXIST=ok) ! PURE functions can't do any I/O
 IF (.NOT. ok) THEN
   PRINT*,'error in geometric_mod: specified file does not exist: ' // TRIM(inputfile)
   STOP
 ENDIF

! Open the file for reading:
 OPEN(UNIT=fid,FILE=TRIM(inputfile),STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ierr)
 IF (ierr/=0) PRINT*, 'problems during openning file :' // TRIM(inputfile)

 ! READ parameters in the order proposed in the input file:
 
!-------------- nv ------------------:
 READ(UNIT=fid, FMT=*,IOSTAT=ierr) pp%nv
 IF(ierr/=0) THEN
   PRINT*,'error during reading number of vertices in geometric_mod.'
   STOP
 ENDIF
!-----!-----!-----!-(xv,yv)----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
 DO i=1,pp%nv
   READ(UNIT=fid, FMT=*, IOSTAT=ierr) pp%xv(i), pp%yv(i)
   IF (ierr/=0) THEN
     PRINT*,'error during reading vertices of the polyprism in geometric_mod.'
     STOP
   ENDIF 
 ENDDO
 !-------------- z1 ------------------:
READ(UNIT=fid, FMT=*,IOSTAT=ierr) pp%z1
 IF (ierr/=0) THEN
   PRINT*,'error during reading top and bottom of the polyprism in geometric_mod.'
   STOP
 ENDIF
 !-------------- z2 ------------------:
READ(UNIT=fid, FMT=*,IOSTAT=ierr) pp%z2
 IF (ierr/=0) THEN
   PRINT*,'error during reading top and bottom of the polyprism in geometric_mod.'
   STOP
 ENDIF
!-------------- ppval ------------------: 
 READ(UNIT=fid, FMT=*,IOSTAT=ierr) pp%ppval
 IF(ierr/=0) THEN
   PRINT*,'error during reading value of physical property in geometric_mod.'
   STOP
 ENDIF
!-------------- props ------------------:
READ(UNIT=fid, FMT=*,IOSTAT=ierr) pp%props
 IF(ierr/=0) THEN
   PRINT*,'error during reading type of physical property in geometric_mod.'
   STOP
 ENDIF
!-------------- inc --------------------:
READ(UNIT=fid, FMT=*,IOSTAT=ierr) pp%inc
 IF(ierr/=0) THEN
   PRINT*,'error during reading inclination in geometric_mod.'
   STOP
 ENDIF
!-------------- dec --------------------:
READ(UNIT=fid, FMT=*,IOSTAT=ierr) pp%dec
 IF(ierr/=0) THEN
   PRINT*,'error during reading declination in geometric_mod.'
   STOP
 ENDIF 
END SUBROUTINE geometric_read_inputs
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
END MODULE geometric_mod
