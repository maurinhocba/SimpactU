    MODULE classNodal_Point

    ! "NodalPoint" extended class type definition
    ! Specialization of "Point"

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE classPoint, only: Point

    IMPLICIT NONE
    PUBLIC

    TYPE, EXTENDS(Point), Public ::  Nodal_Point ! extended type of "Point"
        !Attributes:
        INTEGER            :: Label ! global index of the nodal point
        !INTEGER            :: Set_Position  ! Defines the position in a nodes Set. Exclusive for WAKE use (Replaces the label as for sgmnt-node match)
        DOUBLE PRECISION, DIMENSION(3) :: InducedSpeed  ! Speed induced by grid and wake vorticity at nodal point

    CONTAINS ! methods' names

    PROCEDURE, Public, NoPass :: NPSet_Allocation
    PROCEDURE, Public, NoPass :: NPSet_Generation
    PROCEDURE, Public, NoPass :: NPSet_Printing

    END TYPE Nodal_Point

    CONTAINS ! ===========================================================

    SUBROUTINE NPSet_Allocation (NumberOfNodes, NPoint_Set)  ! Set allocation (Nodal_Point allocatable array)

    !Input Variables
    INTEGER, INTENT(IN)                  :: NumberOfNodes

    !Output Variables
    type(NODAL_POINT), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: NPoint_Set ! Array set variable for allocation

    ALLOCATE (NPoint_Set (NumberOfNodes))

    END SUBROUTINE NPSet_Allocation

    SUBROUTINE  NPSet_Generation (NPSet, NodesArray)
    ! Nodal Points Set generation from NodesArray.

    !Input Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:), INTENT(IN)    :: NodesArray

    !Output variables
    type(NODAL_POINT),  DIMENSION (:), INTENT(INOUT)  :: NPSet

    !Variables Locales
    INTEGER                                           :: i,NumberOfNodes

    NumberOfNodes=size(NPSet)

    do i=1,NumberOfNodes
        NPSet(i)%Label=i
        NPSet(i)%xyz=NodesArray(:,i)
        NPSet(i)%xyz_loc=NodesArray(:,i)
        NpSet(i)%InducedSpeed = [0.0,0.0,0.0]
    end do

    END SUBROUTINE NPSet_Generation

    SUBROUTINE NPSet_Printing (NPoint_Set, FileUnit) ! Prints a Grid's NodalPointSet
    ! File unit specification is needed

    INTEGER, INTENT(IN)                           :: FileUnit
    type(NODAL_POINT),  DIMENSION (:), INTENT(IN) :: NPoint_Set

    INTEGER             :: NumberOfNodes,i

    write (FileUnit,'(A21,/)') 'NodalPoint_Set Data: '
    write (FileUnit,'(A50,/)') 'Label:|   X:   |   Y:   |   Z:  |     InducedSpeed      |'

    NumberOfNodes=size(NPoint_Set)

    DO i=1,NumberOfNodes
        write (FileUnit,60) NPoint_Set(i)%Label, NPoint_Set(i)%xyz, Npoint_Set(i)%InducedSpeed
    end do

60  FORMAT (I5,6F15.4)

    END SUBROUTINE  NPSet_Printing


    END MODULE classNodal_Point