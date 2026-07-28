    MODULE classVortex_Segment

    ! "Vortex_Segments" class definition

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE P_S_N

    IMPLICIT NONE
    PUBLIC


    TYPE, Public :: Vortex_Segment
        !Attributes:
        INTEGER                            :: Grid  ! Grid Owner ID
        INTEGER                            :: Label  ! Global ID
        INTEGER, DIMENSION(2)              :: SegmentNodes ! connectivity list, first-last. Indicates Nodes Labels which may differ from their position in Sets.
        INTEGER                            :: nSharingPanels ! number of panels sharing the segment
        INTEGER, ALLOCATABLE, DIMENSION(:) :: SharingPanels ! panels sharing the segment
        DOUBLE PRECISION                               :: SgmntCirc  ! Circulation of the segment
        DOUBLE PRECISION, DIMENSION (3)                :: Sgmnt_Cond ! (a,b,c) Defines de condition of the segment, edge-convective(b), edge-nonconvective(a), inner, processed(c). For DeltaCP calculation
        DOUBLE PRECISION                               :: SharedArea ! Sum of the Sharing Panels Areas
        DOUBLE PRECISION                               :: Length     ! Length of the segment

    CONTAINS ! methods' names

    PROCEDURE, Public, NoPass :: VSSet_Allocation
    PROCEDURE, Public, NoPass :: VSSet_Generation
    PROCEDURE, Public, NoPass :: VSSet_Printing
    PROCEDURE, Public, NoPass :: VSSet_Length

    END TYPE Vortex_Segment

    CONTAINS ! ==================================================================

    SUBROUTINE VSSet_Allocation (NumberOfSegments, Vtx_Sgmnt_Set)  ! Set allocation (VtxSgmnt allocatable array)

    INTEGER, INTENT(IN)                    :: NumberOfSegments

    type(VORTEX_SEGMENT), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: Vtx_Sgmnt_Set ! Array set variable for allocation

    ALLOCATE (Vtx_Sgmnt_Set (NumberOfSegments))

    Vtx_Sgmnt_Set(:)%SgmntCirc = 0.0
    Vtx_Sgmnt_Set(:)%Label = 0

    END SUBROUTINE VSSet_Allocation

    SUBROUTINE  VSSet_Generation (GridLabel, IENArray, VSSet, NumberOfSegments, ISNArray, IESArray, ShPanels_Links, ShPanels_Storage, PanelType)
    ! Vortex Segment Set generation from Segment_Arrays, works with P_S_N Module


    !INPUT VARIABLES
    INTEGER, INTENT(IN)                               :: GridLabel
    INTEGER, DIMENSION (:,:), INTENT(IN)              :: IENArray
    INTEGER, INTENT(IN)                               :: NumberOfSegments
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(IN) :: ISNArray, IESArray
    INTEGER, ALLOCATABLE, DIMENSION (:), INTENT(IN)   :: ShPanels_Links, ShPanels_Storage
    INTEGER, ALLOCATABLE, DIMENSION (:), INTENT(IN)   :: PanelType

    !OUTPUT VARIABLES
    type(VORTEX_SEGMENT), DIMENSION (:), INTENT(INOUT)  :: VSSet

    !LOCAL VARIABLES
    !INTEGER                               :: NumberOfSegments
    !INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ISNArray, IESArray
    !INTEGER, ALLOCATABLE, DIMENSION (:)   :: ShPanels_Links, ShPanels_Storage
    !INTEGER, ALLOCATABLE, DIMENSION (:)   :: PanelType
    INTEGER  :: i, j, NShPanels

    !call Panel_Type (IENArray, PanelType)
    !call Segments (GridLabel, IENArray, PanelType, NumberOfSegments, ISNArray, IESArray,ShPanels_Links, ShPanels_Storage)


    !NumberOfSegments=size(VSSet)

    do i=1,NumberOfSegments
        VSSET(i)%Grid = GridLabel
        VSSet(i)%Label = i
        VSSet(i)%SegmentNodes = ISNArray(:,i)
        if (i==1) then
            NShPanels = ShPanels_Links(i)
        else
            NShPanels = ShPanels_Links(i)-ShPanels_Links(i-1)
        end if
        VSSet(i)%nSharingPanels=NShPanels

        ALLOCATE ( VSSet(i)%SharingPanels (NShPanels) )
        do j=1,NShPanels
            VSSet(i)%SharingPanels(j) = ShPanels_Storage(ShPanels_Links(i)-j+1)
        end do

        VSSet(i)%Sgmnt_Cond(1) = 1.0
        VSSet(i)%Sgmnt_Cond(2) = 0.0
        VSSet(i)%Sgmnt_Cond(3) = 0.0   ! Sets a default value for further identification = not processed
    end do

    END SUBROUTINE VSSet_Generation

    SUBROUTINE VSSet_Length (VSSet,NPSet)
    ! Calculates the length of every VtxSegment of the set, defined as the 2-norm of the vector joining the nodes of the segment.
    ! VSSet and NPSet must be from the same grid.

    !Input Variables
    type(NODAL_POINT), DIMENSION (:), INTENT(IN)  :: NPSet !Nodal points Set for coordinates extraction

    !Output Variables
    type(VORTEX_SEGMENT), DIMENSION (:), INTENT(INOUT)  :: VSSet

    !Local Variables
    INTEGER :: i !Loop index
    INTEGER :: nVtxSgmnts !Number of Vortex Segments of the set. Extraction variable
    DOUBLE PRECISION, DIMENSION (3) :: N1, N2 !Nodes Coordinates. Extraction variable

    !Number of VtxSgmnts
    nVtxSgmnts = size(VSSet)

    DO i=1,nVtxSgmnts
        N1 = NPSet(VSSet(i)%SegmentNodes(1))%xyz
        N2 = NPSet(VSSet(i)%SegmentNodes(2))%xyz

        VSSet(i)%Length = norm2(N1-N2)
    END DO

    END SUBROUTINE VSSet_Length

    SUBROUTINE VSSet_Printing (VSSet, FileUnit) ! VortexSet printing or file writing (if file specified)
    ! Needs a file specification

    INTEGER, INTENT(IN)      :: FileUnit  !File id number for writing

    type(VORTEX_SEGMENT), DIMENSION(:), INTENT(IN) :: VSSet ! Array set variable for Printing

    INTEGER :: NumberOfSegments,i

    NumberOfSegments=size(VSSet)

    write (FileUnit,'(A25,/)') 'Vortex_Segment_Set Data: '

    write (FileUnit,'(A98,/)') 'Grid: | Label: |   Nodes:   |     Circ:    | NShPanels: |   SharingPanels:   |   Sgmt_Cond: |'

    do i=1,NumberOfSegments
        write (FileUnit, *) VSSet(i)%Grid,VSSet(i)%Label,VSSet(i)%SegmentNodes,VSSet(i)%SgmntCirc,VSSet(i)%nSharingPanels,VSSet(i)%SharingPanels,VSSet(i)%Sgmnt_Cond
    end do

60  FORMAT (4I7,F15.4,I3,'*',2F5.2) ! Increse space for larger grids. ---------------------XXXXXX

    END SUBROUTINE VSSet_Printing



    END MODULE classVortex_Segment
