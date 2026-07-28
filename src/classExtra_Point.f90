MODULE classExtra_Point

    ! "Extra_Point" extended class type definition

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE classPoint, only: Point

    IMPLICIT NONE
    PUBLIC

    TYPE, EXTENDS(Point), Public ::    Extra_Point ! extended type of "Point"
        !Attributes:
        DOUBLE PRECISION, dimension (3) ::    Velocity		  ! DOUBLE PRECISION velocity at Extra point
        DOUBLE PRECISION				::    VelocityMag	  ! DOUBLE PRECISION velocity magnitude

    CONTAINS ! methods' names

    PROCEDURE, Public, NoPass :: ExtraPoint_ReadLoad
    PROCEDURE, Public, NoPass :: ExtPt_Vel

    END TYPE Extra_Point

    CONTAINS ! ==================================================================

SUBROUTINE ExtraPoint_ReadLoad(EP_FileName, NumberOfExtPts, ExtPoint_Set)
    !Reads a Coordinates file for Extra Points and generates, allocates and loads the correspondig set.

    !Input Variables
    CHARACTER(len=20), INTENT(IN)   :: EP_FileName ! Extra Points coordinates file name.
    INTEGER, INTENT(IN)             :: NumberOfExtPts  !Number of Extra Points

    !Output Variables
    type(Extra_Point), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ExtPoint_Set

    !Local Variables
    INTEGER    :: ios !File check variable
    INTEGER    :: i !Loop index
    DOUBLE PRECISION, DIMENSION(3) :: NodeInLine ! File line extraction variable

    !Open File
    open (7, FILE=EP_FileName, IOSTAT=ios, ACTION='read')

    ! Opening status check
    if (ios /= 0) then

        print * , 'File opening error - ExtraPoints File: ',EP_FileName
        pause

    end if

    print * , 'EXTPTS File: ',EP_FileName
    pause
    !File reading - Coordinates
    do i=1,NumberOfExtPts

        read (7,*,IOSTAT=ios) NodeInLine
        print * , 'EXTPTS nil: ', NodeInLine, '/ ', i
        ExtPoint_Set(i)%xyz(1) = NodeInLine(1)
        ExtPoint_Set(i)%xyz(2) = NodeInLine(2)
        ExtPoint_Set(i)%xyz(3) = NodeInLine(3)

        print * , 'EXTPTS xyz: ',  ExtPoint_Set(i)%xyz, '/ '
        ExtPoint_Set(i)%xyz_loc = ExtPoint_Set(i)%xyz

        ExtPoint_Set(i)%Velocity = [0.0, 0.0, 0.0] ! Velocity is initialized in zero.

    end do
    pause
    !Check print
    do i=1,NumberOfExtPts
        print *,ExtPoint_Set(i)%xyz(:),'/'
    end do


    !File Close
    close (UNIT=7)

END SUBROUTINE ExtraPoint_ReadLoad

SUBROUTINE ExtPt_Vel (ExtPoint_Set, V_inf)
    ! Calculates the DOUBLE PRECISION velocity (V_inf + Induced) at the extra point and the magnitude of the velocity vector.

    !Input Variables
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: V_inf

    !Output Variables
    type(Extra_Point), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ExtPoint_Set

    !Local Variables
    INTEGER :: NumberOfExtPts
    INTEGER :: i


    !Size extraction
    NumberOfExtPts = size(ExtPoint_Set)

    DO i=1,NumberOfExtPts

        ExtPoint_Set(i)%Velocity = ExtPoint_Set(i)%Velocity + V_inf  !Adds free srteam velocity
        ExtPoint_Set(i)%VelocityMag = norm2(ExtPoint_Set(i)%Velocity)! Calulates the vector magnitude

    END DO


END SUBROUTINE ExtPt_Vel

END MODULE classExtra_Point