    MODULE classConvection_Line

    ! CONVECTION LINE class definition

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza


    USE classGrid, only: Grid
    USE classVortex_Segment, only: Vortex_Segment
    USE classNodal_Point, only: Nodal_Point
    USE classPanel, only: Panel

    IMPLICIT NONE
    PUBLIC


    TYPE, Public :: Convection_Line

        INTEGER                            :: Owner ! GRID ID. Identifies the grid which holds the line.
        INTEGER                            :: CLineLabel ! ConvecLine ID, sets its reference
        LOGICAl                            :: isClosed ! If true, the CL starts and ends in the same node. Otherwise the line is open
        INTEGER                            :: nNodalPoints ! number of nodal points
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NodalPointSet ! (Set) Nodal points of the ConvecLine
        INTEGER                            :: nVortexSegments ! number of vortex segments
        INTEGER, ALLOCATABLE, DIMENSION(:) :: VortexSegmentSet ! (Set) Vortex segments of the ConvecLine
        CHARACTER                          :: extraData  ! for misc purposes

    CONTAINS ! methods' names

    PROCEDURE, Public, NoPass :: ConvecLine_Reading_VS
    PROCEDURE, Public, NoPass :: ConvecLine_Reading_N
    PROCEDURE, Public, NoPass :: ConvecLine_Loading
    PROCEDURE, Public, NoPass :: ConvecLine_Printing
    PROCEDURE, Public, NoPass :: Segment_Condition

    END TYPE Convection_Line

    CONTAINS ! ==================================================================

    SUBROUTINE ConvecLine_Reading_N (FileName, GRD, CL_VtxSgmnts, nVtxSgmnts)
    !Reads the vtx segments ID set of the Convection Line. Produces an array of integers without further process,
    !also returns the number of Vtx segments.
    ! THE SET UP OF THIS SUBROUTINE DEPENDS STRICTLY OF THE PREPROCESS, DEFAULT IS AS FOLLOWS. Reading from file.

    !Input Variables
    character (len=20),INTENT(IN)                   :: FileName
    type(Grid), INTENT(IN)                          :: GRD

    !Output Varialbes
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: CL_VtxSgmnts
    INTEGER, INTENT(OUT)                            :: nVtxSgmnts

    !Local Variables
    INTEGER :: i,j,k      !Loop index
    INTEGER :: nNodes     !Number of nodes listed in the Convect Line
    INTEGER :: NodeID     !Node ID extraction Variable
    INTEGER :: ios        !Opening status check variable
    INTEGER :: n          !Variable counter
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NodeIDList !Extraction array

    open (UNIT=3, FILE=FileName, IOSTAT=ios, ACTION='read')
    ! Opening status check
    if (ios /= 0) then

        print * , 'File opening error - ConvectionLine File: ', FileName
        pause
    else
        print * ,'ConvectionLine File opening... '
        pause
    end if


    !Nodes Read
    read (3,10) nNodes
    print *, 'nN ', nNodes
    !CL_VtxSgmnts Allocation

    ALLOCATE (CL_VtxSgmnts(nNodes-1))
    ALLOCATE (NodeIDList(nNodes))

    DO i=1,nNodes
        read (3,20) NodeID
        print *, 'nID', NodeID
        NodeIDList(i)=NodeID
    END DO

    n = 0 ! Set counter to zero

    !CL_Nodes Read
    DO j=1,GRD%nVortexSegments
        DO i=1,nNodes
            if (NodeIDList(i) == GRD%VortexSegmentSet(j)%SegmentNodes(1)) then   ! Links Nodes and VtxSgmnts.
                DO k=1,nNodes
                    if (GRD%VortexSegmentSet(j)%SegmentNodes(2) == NodeIDList(k)) then
                        n = n+1
                        CL_VtxSgmnts(n) = GRD%VortexSegmentSet(j)%Label
                        exit
                    end if
                END DO
                exit
            else if (NodeIDList(i) == GRD%VortexSegmentSet(j)%SegmentNodes(2)) then
                DO k=1,nNodes
                    if (GRD%VortexSegmentSet(j)%SegmentNodes(1) == NodeIDList(k)) then
                        n = n+1
                        CL_VtxSgmnts(n) = GRD%VortexSegmentSet(j)%Label
                        exit
                    end if
                END DO
                exit
            end if


        END DO
    END DO

    print *, 'CLVTXS ', CL_VtxSgmnts(:), 'nVtx:    ',n
    nVtxSgmnts = n

    !print*, 'CL_VtxSgmnts', CL_VtxSgmnts
    !   WHEN SEVERAL CONVECTION LINES ARE READ, CROSS REPETITION FOR VORTEX SEGMENTS' ID MUST BE VERIFIED. ------<<<<>>>>>

    close (3)

    deallocate (NodeIDList)

10  FORMAT (//,1i6,/)
20  FORMAT (1i6)
30  FORMAT (/,1i6)

    END SUBROUTINE ConvecLine_Reading_N


    SUBROUTINE ConvecLine_Reading_VS (FileName, CL_VtxSgmnts, nVtxSgmnts)
    !Reads the vtx segments ID set of the Convection Line. Produces an array of integers without further process,
    !also returns the number of Vtx segments.
    ! THE SET UP OF THIS SUBROUTINE DEPENDS STRICTLY OF THE PREPROCESS, DEFAULT IS AS FOLLOWS. Reading from file.

    !Input Variables
    character (len=40),INTENT(IN)                   :: FileName

    !Output Varialbes
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: CL_VtxSgmnts
    INTEGER, INTENT(OUT)                            :: nVtxSgmnts
    !INTEGER, INTENT(OUT)                            :: CLineLabel

    !Local Variables
    INTEGER :: i,j        !Loop index
    INTEGER :: VtxS     !Extraction Variable
    INTEGER :: ios      !Opening status check variable


    open (UNIT=3, FILE=FileName, IOSTAT=ios, ACTION='read')
    ! Opening status check
    if (ios /= 0) then

        print * , 'File opening error'
        pause

    end if

    !nVtxSgmnts Read
    !read (3,30) CLineLabel    !READS CL LABEL FROM FILE ---------------------- remove!!!!!!
    read (3,10) nVtxSgmnts

    !CL_VtxSgmnts Allocation

    ALLOCATE(CL_VtxSgmnts(nVtxSgmnts))

    !CL_VtxSgmnts Read
    DO i=1,nVtxSgmnts
        read (3,20) VtxS
        CL_VtxSgmnts(i) = VtxS
        DO j=1,i-1

            if (VtxS == CL_VtxSgmnts(j)) then   ! Checks for repeted vtx segments in the list.
                error stop "ERROR: Segments in current Convection Line are repeted."   !--------->>>> CREATE ERROR OUTPUT FILE !!!!!
            end if
        END DO

    END DO
    !print*, 'CL_VtxSgmnts', CL_VtxSgmnts
    !   WHEN SEVERAL CONVECTION LINES ARE READ, CROSS REPETITION FOR VORTEX SEGMENTS' ID MUST BE VERIFIED. ------<<<<>>>>>

    close (3)


10  FORMAT (/,1i6,/)
20  FORMAT (1i6)
30  FORMAT (/,1i6)

    END SUBROUTINE ConvecLine_Reading_VS

    SUBROUTINE ConvecLine_Loading (GRD, CL_VtxSgmnts, nVtxSgmnts, CLineLabel, ConvecLine)
    !Loads and processes data read from file to Convection_Line class.

    ! Input Variables
    type(GRID), INTENT (IN)                        :: GRD
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: CL_VtxSgmnts
    INTEGER, INTENT(IN)                            :: nVtxSgmnts
    INTEGER, INTENT(IN)                            :: CLineLabel

    !Output Variables
    type(Convection_Line), INTENT(OUT)    :: ConvecLine

    !Local Variables
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CL_Nodes   ! Primary nodes extraction array
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: FrontNodes ! Temp Array for nodes sort
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: BackNodes  ! Temp Array for nodes sort
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: FrontVS    ! Temp Array for vtx sgmnts sort
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: BackVS     ! Temp Array for vtx sgmnts sort
    INTEGER                              :: i,j        ! Loop Indexes
    INTEGER                              :: FrntCNT    ! Aux Counter
    INTEGER                              :: FinalNode, InitialNode !First and Last nodes of the CL

    ALLOCATE (FrontNodes(nVtxSgmnts))
    ALLOCATE (FrontVs (nVtxSgmnts))
    ALLOCATE (CL_Nodes(2,nVtxSgmnts))

    !Nodes primary extraction
    DO i=1,nVtxSgmnts
        CL_Nodes(1,i) = GRD%VortexSegmentSet(CL_VtxSgmnts(i))%SegmentNodes(1)
        CL_Nodes(2,i) = GRD%VortexSegmentSet(CL_VtxSgmnts(i))%SegmentNodes(2)
    END DO
    !print*,'CL_Nodes',CL_Nodes

    !Nodes and VS sorting
    FrntCNT = 0
    FinalNode = 0
    FrontNodes(1) = CL_Nodes (2,1)
    FrontVS(1) = CL_VtxSgmnts (1)
    DO i = 2,nVtxSgmnts
        if (FinalNode == 0) then
            DO j = 2,nVtxSgmnts
                if (CL_Nodes(1,j)==FrontNodes(i-1)) then !.and. CL_Nodes(2,j)/=0) then
                    FrontVS(i) = CL_VtxSgmnts(j)
                    FrontNodes(i) = CL_Nodes(2,j)
                    CL_Nodes(2,j) = 0
                    FrntCNT = FrntCNT + 1
                    exit
                else if (CL_Nodes(2,j)==FrontNodes(i-1)) then !.and. CL_Nodes(1,j)/=0) then
                    FrontVS(i) = -CL_VtxSgmnts(j)
                    FrontNodes(i) = CL_Nodes(1,j)
                    CL_Nodes(1,j) = 0
                    FrntCNT = FrntCNT + 1
                    exit
                else if (j == nVtxSgmnts) then
                    FinalNode = FrontNodes(i-1)
                end if
            END DO
        end if
    END DO

    !print*,'FrntCNT ',FrntCNT
    !print*,'FinalNode ',FinalNode

    if (FinalNode==0)then
        FinalNode = FrontNodes(nVtxSgmnts)
    end if

    ALLOCATE (BackNodes(nVtxSgmnts-FrntCNT))
    ALLOCATE (BackVs (nVtxSgmnts-FrntCNT-1))

    BackNodes(1) = CL_Nodes (1,1)
    InitialNode = 0
    DO i = 2,nVtxSgmnts-FrntCNT
        if (InitialNode == 0) then
            DO j = 2,nVtxSgmnts
                if (CL_Nodes(1,j)==BackNodes(i-1)) then!  .and. CL_Nodes(2,j)/=0) then
                    BackVS(i-1) = -CL_VtxSgmnts(j)
                    BackNodes(i) = CL_Nodes(2,j)
                    CL_Nodes(2,j) = 0
                    exit
                else if (CL_Nodes(2,j)==BackNodes(i-1)) then!  .and. CL_Nodes(1,j)/=0) then
                    BackVS(i-1) = CL_VtxSgmnts(j)
                    BackNodes(i) = CL_Nodes(1,j)
                    CL_Nodes(1,j) = 0
                    exit
                else if (j == nVtxSgmnts)then
                    InitialNode = BackNodes(i-1)
                end if
            END DO
        end if
    END DO
    print*,'InitialNode ',InitialNode
    if (InitialNode==0) then
        InitialNode = BackNodes(nVtxSgmnts-FrntCNT)
    end if

    !print*,'FrontVS ',FrontVS
    !print*,'BackVS ',BackVS
    !print*,'FrontN ',FrontNodes
    !print*,'BackN ',BackNodes

    !CL continuity verification, execution must stop if false.
    if (BackNodes(nVtxSgmnts-FrntCNT) == 0) then    !Assert condition.
        error stop "ERROR: Convection line is not continuous."   !--------->>>> CREATE ERROR OUTPUT FILE !!!!!
    end if

    ! Convection Line clasification
    if (InitialNode==FinalNode) then
        ConvecLine%isClosed = .true.   !CL is CLOSED
        ConvecLine%nNodalPoints = nVtxSgmnts
    else
        ConvecLine%isClosed = .false.  !CL is OPEN
        ConvecLine%nNodalPoints = nVtxSgmnts + 1
    end if

    !Final Load
    print*,'FinalLoad NPS: '
    ALLOCATE (ConvecLine%NodalPointSet (ConvecLine%nNodalPoints))
    ALLOCATE (ConvecLine%VortexSegmentSet (nVtxSgmnts))

    DO i=1,FrntCNT+1
        ConvecLine%NodalPointSet(nVtxSgmnts-FrntCNT+i) = FrontNodes(i)
        ConvecLine%VortexSegmentSet(nVtxSgmnts-FrntCNT+i-1) = FrontVs(i)
    END DO

    if (ConvecLine%isClosed == .false.) then
        DO i=1,(nVtxSgmnts-FrntCNT)
            ConvecLine%NodalPointSet(i) = BackNodes(nVtxSgmnts-FrntCNT-i+1)
        END DO
        DO i=1,(nVtxSgmnts-FrntCNT-1)
            ConvecLine%VortexSegmentSet(i) = BackVs(nVtxSgmnts-FrntCNT-i)
        END DO
    end if

    ConvecLine%Owner = GRD%GridLabel
    ConvecLine%CLineLabel = CLineLabel  ! --------> NUMERATION SYSTEM
    ConvecLine%nVortexSegments = nVtxSgmnts

    print*,'FinalLoad---- '
    print*,'ConvLine Segments: ', ConvecLine%VortexSegmentSet
    pause

    END SUBROUTINE ConvecLine_Loading

    SUBROUTINE ConvecLine_Printing (ConvecLine, FileUnit) ! VortexSet printing or file writing (if file specified)
    ! Needs a file specification

    !Input Variables
    INTEGER, INTENT(IN)      :: FileUnit  !File id number for writing
    type(CONVECTION_LINE),INTENT(IN) :: ConvecLine ! Convection Line variable for Printing

    !Local Variables
    INTEGER :: i   !Loop index

    write (FileUnit,'(A23,/)') 'Convection_Line Data: '

    write (FileUnit,10) 'Label: ', ConvecLine%CLineLabel
    !write (FileUnit,11) 'Extra Data: ', ConvecLine%extraData  ----> Extra Data is not defined---
    write (FileUnit,12) 'Is Closed?: ', ConvecLine%isClosed
    write (FileUnit,13) 'Number Of Nodes: ', ConvecLine%nNodalPoints
    write (FileUnit,14) 'Number Of Segments: ', ConvecLine%nVortexSegments

    write (FileUnit,'(A7)') 'Nodes:'
    do i=1,ConvecLine%nNodalPoints
        write (FileUnit, 15) ConvecLine%NodalPointSet(i)
    end do
    write (FileUnit,'(A15)') 'VortexSegments:'
    do i=1,ConvecLine%nVortexSegments
        write (FileUnit, 16) ConvecLine%VortexSegmentSet(i)
    end do

    write (FileUnit,'(/)')

    !FORMATS
10  FORMAT (A6,I7)
    !11 FORMAT (A12,A15) ----> Extra Data is not defined---
12  FORMAT (A12,L1)
13  FORMAT (A17,I9)
14  FORMAT (A20,I9)
15  FORMAT (I5) ! Increse space for larger grids. ---------------------------------XXXXXX
16  FORMAT (I5) ! Increse space for larger grids. ---------------------------------XXXXXX

    END SUBROUTINE ConvecLine_Printing

    SUBROUTINE Segment_Condition (GRD, ConvLine)
    !Determines the condition of a grid Vortex Segment respecting the "edge condition",
    !i.e. "Inner Segment", "Convective Edge", "Non-Convective Edge".
    !Used in DeltaCP calculation

    !Input Variables
    type(CONVECTION_LINE), INTENT(IN) :: ConvLine

    !Output Variables
    type(GRID), INTENT(INOUT)         :: GRD

    !Local Variables
    INTEGER :: i, j !Loop indexes
    !INTEGER :: GrdVSLabel, ClnVSLabel ! Labels extraction variables

    DO i=1,GRD%nVortexSegments
        IF (GRD%VortexSegmentSet(i)%Sgmnt_Cond(3) == 0.0) THEN
            !GRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 1.0
            !GRD%VortexSegmentSet(i)%Sgmnt_Cond(2) = 0.0

            !Edge Identification
            if (GRD%VortexSegmentSet(i)%nSharingPanels == 1) then
                GRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 2.0

                !Convective Identification
                DO j=1,ConvLine%nVortexSegments

                    if (GRD%VortexSegmentSet(i)%Label == abs(ConvLine%VortexSegmentSet(j))) then
                        GRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 1.0
                        GRD%VortexSegmentSet(i)%Sgmnt_Cond(2) = 1.0
                        GRD%VortexSegmentSet(i)%Sgmnt_Cond(3) = 1.0
                        exit
                        !else if (j == ConvLine%nVortexSegments) then
                        !    GRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 2.0
                    end if
                END DO
            end if
        END IF

    END DO

    END SUBROUTINE Segment_Condition

    END MODULE classConvection_Line