    MODULE classWake

    ! Wake abstract data type definition

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE classVortex_Segment
    USE classNodal_Point
    !USE classWake_Panel
    USE classGrid
    USE classConvection_Line
    !USE P_S_N

    IMPLICIT NONE
    PUBLIC


    TYPE, Public :: Wake

        CHARACTER(len=20)                               :: WakeName
        INTEGER                                         :: WakeLabel ! Wake ID. Sets the reference for components (nodes, panels, vtxsgments...)
        CHARACTER (len=15)                              :: ExtraData  ! for miscellaneous purposes
        INTEGER                                         :: Owner ! Grid ID from which the wake is shed
        INTEGER                                         :: WakeSize ! Defines the lenght of the Wake in Nodes/Segments Block numbers.
        INTEGER                                         :: LLStride ! Set the position of the last loaded stride in the wake. If wake is full *=Posic, else *=-1
        INTEGER                                         :: NodesStride ! Number of nodes of a stride, size of the nodes block
        INTEGER                                         :: SgmntsStride ! Number of segments of a stride, size of the segments block
        type(Convection_Line)                           :: ConvecLine ! Convection Line from where the wake is shed
        !INTEGER                                         :: ConvecLine ! Convection Line ID from where the wake is shed
        INTEGER                                         :: nNodalPoints ! number of nodal points
        type(NODAL_POINT), ALLOCATABLE, DIMENSION(:)    :: NodalPointSet ! (Set) Nodal points of the Wake
        INTEGER                                         :: nVortexSegments ! number of vortex segments
        type(VORTEX_SEGMENT), ALLOCATABLE, DIMENSION(:) :: VortexSegmentSet ! (Set) Vortex segments of the Wake
        !INTEGER                                         :: nWkPanels ! number of panels
        !type(WAKE_PANEL), ALLOCATABLE, DIMENSION(:)     :: WkPanelSet  ! (Set) Panels of the Wake

    CONTAINS ! methods' names

    !PROCEDURE, Public, Nopass :: IndSpeed_atNode
    !PROCEDURE, Public, Nopass :: Wake_NodesIndSpeed
    PROCEDURE, Public, Nopass :: Wake_Allocation
    PROCEDURE, Public, Nopass :: Wake_Initialization
    PROCEDURE, Public, Nopass :: Wake_NewStep
    PROCEDURE, Public, Nopass :: Wake_Printing
    !PROCEDURE, Public, Nopass :: CONVECTION
    PROCEDURE, Public, Nopass :: Wake_Influence_IndSpeed
    PROCEDURE, Public, Nopass :: Wake_Reading

    END TYPE Wake

    CONTAINS ! ==================================================================

    SUBROUTINE Wake_Printing (WKE, FileUnit)
    !Prints a Grid into a file. Requires the file to be opened and available for writing.

    !Input Variables
    INTEGER, INTENT(IN)    :: FileUnit
    type(WAKE), INTENT(IN) :: WKE

    !Local Variables


    !Grid General Data Printing

    write (FileUnit,'(A19,/)') 'WAKE General Data: '
    write (FileUnit,10) 'Name: ', WKE%WakeName
    write (FileUnit,11) 'Label: ', WKE%WakeLabel
    !write (FileUnit,12) 'Extra Data: ', WKE%ExtraData  ----> Extra Data is not defined---
    write (FileUnit,13) 'Owner: ', WKE%Owner
    write (FileUnit,14) 'Number Of Nodes: ', WKE%nNodalPoints
    write (FileUnit,15) 'Number Of Segments: ', WKE%nVortexSegments
    write (FileUnit,'(/)')

    ! Nodal Point Set Printing
    call NPSet_Printing (WKE%NodalPointSet, FileUnit)

    ! Vortex Segments Set Printing
    call VSSet_Printing (WKE%VortexSegmentSet, FileUnit)

    !FORMATS
10  FORMAT (A6,A)
11  FORMAT (A7,I10)
    !12 FORMAT (A12,A15) ----> Extra Data is not defined---
13  FORMAT (A7,I5)
14  FORMAT (A17,I9)
15  FORMAT (A17,I9)



    END SUBROUTINE Wake_Printing


    !SUBROUTINE CONVECTION (BodySet, nBodies, TStep, StepCnt, V_inf)
    !    ! Defines and organices the convection process for a Grid with ConvectLine into a Wake and goes through every wake in the BodySet.
    !    ! WkeSize and TStep are the number of convection steps that will be stored and their time step, respectivly.
    !    ! Includes the specific routines in this module, asuming that the wake is already allocatad and initiallizated.
    !    ! Assumes that the Wake and its Owners (ConvectLine and Grid) are from the same body.
    !    ! Requires previous Induced Speed at Wake Nodes calculation.
    !
    !!Input Variables
    !INTEGER, INTENT(IN)               :: nBodies
    !INTEGER, INTENT(IN)               :: StepCnt
    !DOUBLE PRECISION, INTENT(IN)                  :: TStep
    !DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: V_inf
    !
    !!Output Variables
    !type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    !
    !!Local Variables
    !type(Nodal_Point), ALLOCATABLE, DIMENSION(:)    :: NewNodes ! Group of nodes generated in each step of convection, forming a "parallel" line from convection line.
    !type(Vortex_Segment), ALLOCATABLE, DIMENSION(:) :: NewSegments !Group of segments genereated according to the 'NewNodes' array.
    !INTEGER     :: i,j,k  !Loop Index
    !
    !
    !DO i=1,nBodies
    !    DO j,1=BodySet(i)%nWakes
    !
    !    WKEj: associate (WKE => BodySet(i)%Wakes(j))
    !
    !        !Owner Search
    !        DO k=1, BDY%nGrids
    !            if (BDY%Grids(k)%GridLabel == WKE%Owner) then
    !                OwAux = k
    !                Ownerj: associate (OWNR => BodySet(i)%Grids(k))
    !                exit
    !            end if
    !        END DO
    !
    !    call Wake_NewStep (WKE%ConvectionLine, OWNR, TStep, NewNodes, NewSegments, V_inf)   ! Generates generic new nodes and segments as arrays
    !
    !    call Wake_Actualization (WKE, NewNodes, NewSegments, StepCnt, WKE%WakeSize, TStep, V_inf) ! Inserts the new nodes and segments in the existing wake and actualices the position and indexes of the previous ones.
    !
    !        end associate WKEj
    !    end associate Ownerj
    ! END DO
    ! END DO
    !
    !END SUBROUTINE CONVECTION

    SUBROUTINE Wake_Influence_IndSpeed (WKE, NodeXyz, WkeIndSpeed)
    !Calculates the induced speed by the wake at a node.
    !Contains the conversion from nodes label to array position for segment connectivities.
    interface
    FUNCTION VtxSgmnt_IndSpeed (NP1, NP2, CP, Gamma)
    !Calculates the speed induced by a vortex segment (from NP1 to NP2) of intensity Gamma, in the control point CP.

    IMPLICIT NONE

    !Input Variables
    DOUBLE PRECISION, DIMENSION (3)  :: NP1, NP2    ! Nodal Points of the Segment, first-last.
    DOUBLE PRECISION, DIMENSION (3)  :: CP          ! Control Point to calculate induced speed in.
    DOUBLE PRECISION                 :: Gamma       ! Circulation of the segment.

    !Output Variables
    DOUBLE PRECISION, DIMENSION (3)  :: VtxSgmnt_IndSpeed   ! Induced speed in the control point.

    !Local Variables
    DOUBLE PRECISION, DIMENSION (3)  :: R1, R2          ! Pointing vectors from CP to NPi.
    DOUBLE PRECISION, PARAMETER      :: PI=3.1415926535 ! Pi.
    DOUBLE PRECISION                 :: L               ! Length of the vortex segment.
    DOUBLE PRECISION                 :: delta = 0.05    ! "Cutoff" radius parameter (VAN GARREL)
    END FUNCTION VtxSgmnt_IndSpeed
    end interface

    !Input Variables
    type(WAKE), INTENT(IN)   :: WKE     !Input wake, for vtxsgmnt set influence
    DOUBLE PRECISION, DIMENSION(3)       :: NodeXyz ! Node where the induced speed is calculated

    !OutputVariables
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: WkeIndSpeed  ! Cumulative Induced Speed

    !LocalVariables
    INTEGER             :: i,j  !Loop Indexes
    INTEGER             :: N1,N2   ! Temp nodes id (Label=Position)
    INTEGER             :: L1, L2  ! Temp nodes lables
    DOUBLE PRECISION, DIMENSION (3) :: xyz1, xyz2 ! Temp nodes coordinates
    INTEGER             :: NodesStride !Size of the block of nodes of the wake (=InitNodes)
    INTEGER             :: CNT

    WkeIndSpeed =  [0.0,0.0,0.0]

    NodesStride = WKE%NodesStride

    ! Wake Vortex Segments
    WkeVS: associate (nWkeVS => WKE%nVortexSegments)
        CNT = 0

        if (WKE%LLStride < 0 .or. WKE%LLStride==WKE%WakeSize-1) then  ! WAKE is not full yet .or. the last stride was loaded at te last position.
            DO j=1,WKE%SgmntsStride*abs(WKE%LLStride+1)+(WKE%SgmntsStride-1)/2
                !DO j=1,nWkeVS
                L1 = WKE%VortexSegmentSet(j)%SegmentNodes(1)
                L2 = WKE%VortexSegmentSet(j)%SegmentNodes(2)

                !print*,'LLStride', L1, L2

                !if ((L1 /= 0) .and. (L2 /= 0)) then
                xyz1 = WKE%NodalPointSet(L1)%xyz
                xyz2 = WKE%NodalPointSet(L2)%xyz
                ! else
                !  exit
                ! end if

                WkeIndSpeed = WkeIndSpeed + VtxSgmnt_IndSpeed (xyz1, xyz2, NodeXyz, WKE%VortexSegmentSet(j)%SgmntCirc)

            END DO

        else
            DO j=1,nWkeVS
                L1 = WKE%VortexSegmentSet(j)%SegmentNodes(1)
                L2 = WKE%VortexSegmentSet(j)%SegmentNodes(2)

                if (L1 <= NodesStride) then
                    CNT = CNT +1
                    N1 = L1
                elseif (L1 <= (WKE%nNodalPoints-((WKE%LLStride+1)*NodesStride))) then
                    CNT = CNT +1
                    N1 = L1 + (WKE%LLStride+1)*NodesStride
                else
                    CNT = CNT +1
                    N1 = L1 - (WKE%WakeSize-(WKE%LLStride+1))*NodesStride
                end if

                if (L2 <= NodesStride) then
                    N2 = L2
                elseif (L2 <= (WKE%nNodalPoints-((WKE%LLStride+1)*NodesStride))) then
                    N2 = L2 + (WKE%LLStride+1)*NodesStride
                else
                    N2 = L2 - (WKE%WakeSize-(WKE%LLStride+1))*NodesStride
                end if

                xyz1 = WKE%NodalPointSet(N1)%xyz
                xyz2 = WKE%NodalPointSet(N2)%xyz

                WkeIndSpeed = WkeIndSpeed + VtxSgmnt_IndSpeed (xyz1, xyz2, NodeXyz, WKE%VortexSegmentSet(j)%SgmntCirc)

            END DO
        end if
    end associate WkeVS

    END SUBROUTINE Wake_Influence_IndSpeed


    !SUBROUTINE Wake_NodesIndSpeed (GRD, WKE)
    !    !Calculates the speed induced by the wake and grid influences at every node of the wake.
    !    !Prior step for convection.
    !
    !!Input Variables
    !type(Grid), INTENT(IN)  :: GRD
    !
    !!Output Variables
    !type(Wake), INTENT(INOUT) :: WKE
    !
    !!Loca Variables
    !INTEGER :: nNodes
    !INTEGER :: i
    !
    !nNodes = WKE%nNodalPoints
    !
    !!DO i=1,nNodes
    !if (WKE%LLStride < 0 .or. WKE%LLStride==WKE%WakeSize-1) then
    !
    !    DO i=1,WKE%NodesStride*abs(WKE%LLStride+1)+WKE%NodesStride
    !        !if (WKE%NodalPointSet(i)%Label == 0) then
    !        !exit
    !        !else
    !        call IndSpeed_atNode (GRD, WKE, WKE%NodalPointSet(i))
    !    !end if
    !    END DO
    !else
    !    DO i=1,nNodes
    !        call IndSpeed_atNode (GRD, WKE, WKE%NodalPointSet(i))
    !    END DO
    !end if
    !
    !
    !END SUBROUTINE Wake_NodesIndSpeed

    SUBROUTINE Wake_Allocation (WKE, ConvectionLine, WkeSize)
    !Allocates memory for wake. It's defined by the number of nodes of the convection line times the number of steps
    !of convection. Information beyond that point is negligible and not considered in further calculations.

    !Input Variables
    type(Convection_Line), INTENT(IN)  :: ConvectionLine
    INTEGER, INTENT(IN)                :: WkeSize

    !Output Variables
    type(Wake), INTENT(INOUT)       :: WKE

    !Local Variables
    INTEGER      :: nCL_nodes
    INTEGER      :: nCL_sgmnts

    !Convection Line and Grid reference

    WKE%ConvecLine = ConvectionLine
    WKE%Owner = ConvectionLine%Owner

    !Segment and nodes count
    nCL_nodes = ConvectionLine%nNodalPoints
    nCL_sgmnts = ConvectionLine%nVortexSegments

    !Wake size definition
    WKE%NodesStride = nCL_nodes
    WKE%SgmntsStride = (nCL_sgmnts*2+1)
    WKE%WakeSize = WkeSize
    Wke%nNodalPoints = nCL_nodes * WkeSize + nCL_nodes
    Wke%nVortexSegments = (nCL_sgmnts*2+1) * WkeSize + nCL_sgmnts   ! The wake must be closed by the last line of segments.

    ! Sets allocation
    !print *, 'VSSet Alloc...', size(Wke%VortexSegmentSet),allocated(Wke%VortexSegmentSet)
    !print *, 'NPSet size...',size(Wke%NodalPointSet)
    !if (allocated(Wke%VortexSegmentSet)) then
    !    print *, 'VSSet DEAlloc...'
    !deallocate(Wke%VortexSegmentSet)
    !end if

    call VSSet_Allocation (Wke%nVortexSegments, Wke%VortexSegmentSet)
    !print *, 'NPSet SizeAlloc...',size(Wke%NodalPointSet),allocated(Wke%NodalPointSet)
    !if (allocated(Wke%NodalPointSet)) then
    !deallocate(Wke%NodalPointSet)
    !end if
    !print *, 'NPSet Alloc...'

    call NPSet_Allocation (Wke%nNodalPoints, Wke%NodalPointSet)
    !print *, 'Alloc END...'
    END SUBROUTINE Wake_Allocation

    SUBROUTINE Wake_Initialization (WKE, ConvectionLine, GRD)
    ! Initialize the wake with the nodes and segments of the convection line. First line of nodes and segments is created
    ! for further convection.

    !Input Variables
    type(Convection_Line), INTENT(IN)  :: ConvectionLine
    type(Grid), INTENT(IN)             :: GRD

    !Output Variables
    type(Wake), INTENT(INOUT)     :: WKE

    !Local Variables
    INTEGER    :: NodeID
    INTEGER    :: i,j  ! Loop indexes

    !Default Posic
    WKE%LLStride = -1  ! Sets the configuration when there's room in the wake for further convection.

    DO i=1,ConvectionLine%nNodalPoints !Copies convection nodes geometric data for first line of wake.
        !print *, 'i 1',i

        NodeID = ConvectionLine%NodalPointSet(i)
        !print *, 'i 2,size',i, size(WKE%NodalPointSet)

        WKE%NodalPointSet(i)%Label = i
        !print *, 'i 3',i

        WKE%NodalPointSet(i)%xyz = GRD%NodalPointSet(NodeID)%xyz
        WKE%NodalPointSet(i)%InducedSpeed = [0.0,0.0,0.0]
    END DO

    !DO i=ConvectionLine%nNodalPoints+1, WKE%nNodalPoints
    !    NodeID = 0
    !    WKE%NodalPointSet(i)%Label = 0
    !    WKE%NodalPointSet(i)%xyz = [0.0,0.0,0.0]
    !    WKE%NodalPointSet(i)%InducedSpeed = [0.0,0.0,0.0]
    !END DO

    DO j=1,ConvectionLine%nNodalPoints-1 !Copies convection segments geometric data for first line of wake.
        WKE%VortexSegmentSet(j)%Grid = WKE%WakeLabel    !Wake ID
        WKE%VortexSegmentSet(j)%Label = j
        WKE%VortexSegmentSet(j)%SegmentNodes(1)=j
        WKE%VortexSegmentSet(j)%SegmentNodes(2)=j+1
        WKE%VortexSegmentSet(j)%SgmntCirc = 0.0  ! Segments are initialized with no circulation, awaiting for later convection
    END DO
    print *, 'INIT END...'
    END SUBROUTINE

    SUBROUTINE Wake_NewStep (ConvectionLine, WkeLabel, GRD, TStep, NewNodes, NewSegments, V_inf)
    ! Generates a step of convection of a wake, from the ConvectionLine of the grid (GRD) acording to the time step (TStep)
    ! specified. Produces a new set of nodes and the corresponding segments as separate structures.

    !Input Variables
    type(Convection_Line), INTENT(IN)  :: ConvectionLine
    type(Grid), INTENT(IN)             :: GRD
    INTEGER                            :: WkeLabel
    DOUBLE PRECISION, INTENT (IN)                  :: TStep
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)     :: V_inf

    !Output Variables
    type(Nodal_Point), ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: NewNodes ! Group of nodes generated in each step of convection, forming a "parallel" line from convection line.
    type(Vortex_Segment), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: NewSegments !Group of segments genereated according to the 'NewNodes' array.

    !Local Variables
    INTEGER              :: i,j !Loop indexes
    DOUBLE PRECISION, DIMENSION(3)   :: IndSpeed ! Induced speed at node. Temp variable
    INTEGER              :: NodeID ! Node label. Temp variable
    INTEGER              :: SgmntID1, SgmntID2 ! VtxSgmnt label. Temp variable
    DOUBLE PRECISION, DIMENSION(3)   :: NewCoords ! New node coordinates.
    INTEGER              :: nVtxSgmnts
    DOUBLE PRECISION                 :: AbsCirc1, AbsCirc2 !Extracts the circulation of a vortex segment with the sign given by Grid
    !DOUBLE PRECISION                 :: beta ! Convection Coeficient, for different sized panels.

    !First step geometry definition: Nodes
    allocate (NewNodes(ConvectionLine%nNodalPoints))

    DO i=1,ConvectionLine%nNodalPoints
        NodeID = ConvectionLine%NodalPointSet(i)
        IndSpeed = GRD%NodalPointSet(NodeID)%InducedSpeed
        NewCoords =  GRD%NodalPointSet(NodeID)%xyz
        NewNodes(i)%Label = i
        NewNodes(i)%xyz = NewCoords
        NewNodes(i)%InducedSpeed = [0.0,0.0,0.0]
    END DO

    !First step GEOMETRY definition: Vtx Segments
    nVtxSgmnts = 2*ConvectionLine%nVortexSegments+1         ! For the Convection Line is nVtxSgmnts = nNodalPoints + 1
    allocate (NewSegments(nVtxSgmnts))

    DO j=1,ConvectionLine%nVortexSegments+1
        NewSegments(j)%Grid = WkeLabel    !Wake ID
        NewSegments(j)%Label = j
        NewSegments(j)%SegmentNodes(1)=j
        NewSegments(j)%SegmentNodes(2)=j-ConvectionLine%nVortexSegments-1
        NewSegments(j)%SgmntCirc = 0.0 ! Segment circulation is set to zero to avoid acumulation mistakes. It is overwritten later.
    END DO

    DO j=ConvectionLine%nVortexSegments+2,nVtxSgmnts
        NewSegments(j)%Grid = WkeLabel
        NewSegments(j)%Label = j
        NewSegments(j)%SegmentNodes(1)=j-ConvectionLine%nVortexSegments-1
        NewSegments(j)%SegmentNodes(2)=j-ConvectionLine%nVortexSegments
        NewSegments(j)%SgmntCirc = 0.0 ! Segment circulation is set to zero to avoid acumulation mistakes. It is overwritten later.
    END DO


    !First step CIRCULATION definition: Vtx Segments
    !beta = (GRD%VortexSegmentSet(abs(ConvectionLine%VortexSegmentSet(1)))%Length)/sqrt(GRD%VortexSegmentSet(abs(ConvectionLine%VortexSegmentSet(1)))%SharedArea)
    NewSegments(1)%SgmntCirc = (GRD%VortexSegmentSet(abs(ConvectionLine%VortexSegmentSet(1)))%SgmntCirc * sign(1,ConvectionLine%VortexSegmentSet(1))) !* beta
    DO j=2,ConvectionLine%nVortexSegments
        SgmntID2 = ConvectionLine%VortexSegmentSet(j)   !Has the sign of the circulation direction
        SgmntID1 = ConvectionLine%VortexSegmentSet(j-1) !Has the sign of the circulation direction
        !beta = (GRD%VortexSegmentSet(abs(SgmntID2))%Length)/sqrt(GRD%VortexSegmentSet(abs(SgmntID2))%SharedArea) ! Calculates the segment convection coefficient L^2/Area
        absCirc2 = GRD%VortexSegmentSet(abs(SgmntID2))%SgmntCirc! * beta !Extracts the circulation of the segment with its original sign, from grid
        !beta = (GRD%VortexSegmentSet(abs(SgmntID1))%Length)/sqrt(GRD%VortexSegmentSet(abs(SgmntID1))%SharedArea) ! Calculates the segment convection coefficient L^2/Area
        absCirc1 = GRD%VortexSegmentSet(abs(SgmntID1))%SgmntCirc! * beta !Extracts the circulation of the segment with its original sign, from grid
        NewSegments(j)%SgmntCirc = (absCirc2 * sign(1,SgmntID2) - absCirc1 * sign(1,SgmntID1))!*1.2 ! Includes the sign defined by the direction of the convection line
    END DO
    !beta = (GRD%VortexSegmentSet(abs(ConvectionLine%VortexSegmentSet(ConvectionLine%nVortexSegments)))%Length)/sqrt(GRD%VortexSegmentSet(abs(ConvectionLine%VortexSegmentSet(ConvectionLine%nVortexSegments)))%SharedArea)
    NewSegments(ConvectionLine%nVortexSegments+1)%SgmntCirc = (- GRD%VortexSegmentSet(abs(ConvectionLine%VortexSegmentSet(ConvectionLine%nVortexSegments)))%SgmntCirc &
        * sign(1,ConvectionLine%VortexSegmentSet(ConvectionLine%nVortexSegments)))! *beta

    !NewSegments(ConvectionLine%nVortexSegments+1)%SgmntCirc = NewSegments(ConvectionLine%nVortexSegments+1)%SgmntCirc &
    !    + GRD%VortexSegmentSet(abs(ConvectionLine%VortexSegmentSet(1)))%SgmntCirc * sign(1,ConvectionLine%VortexSegmentSet(1))

    DO j=ConvectionLine%nVortexSegments+2,nVtxSgmnts
        SgmntID1 = ConvectionLine%VortexSegmentSet(j-ConvectionLine%nVortexSegments-1) !Has the sign of the circulation direction
        !beta = (GRD%VortexSegmentSet(abs(SgmntID1))%Length)/sqrt(GRD%VortexSegmentSet(abs(SgmntID1))%SharedArea) ! Calculates the segment convection coefficient L^2/Area
        absCirc1 = GRD%VortexSegmentSet(abs(SgmntID1))%SgmntCirc !Extracts the circulation of the segment with its original sign, from grid.
        NewSegments(j)%SgmntCirc = (- absCirc1 * sign(1,SgmntID1))! / beta ! Includes the sign defined by the direction of the convection line
    END DO


    END SUBROUTINE Wake_NewStep

    SUBROUTINE Wake_Actualization (WKE, NewNodes, NewSegments, StepCnt, WkeSize, TStep, V_inf)
    ! Inserts the new set of nodes and segments into the existing (already initialized) wake, and actualizes the position of
    ! the existing nodes according to the step characteristics.
    ! ASUMES THAT INDUCED SPEED IS ALREADY CALCULATED AND STORED AT EACH EXISTING NODE OF THE WAKE PRIOR TO EVERY STEP.

    !Input Variables
    type(Nodal_Point), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)    :: NewNodes     ! Set of new nodes generated by 'Wake_NewStep' for insertion and deallocation.
    type(Vortex_Segment), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: NewSegments  ! Set of new segments generated by 'Wake_NewStep' for insertion and deallocation.
    INTEGER, INTENT(IN)                                            :: StepCnt      ! Number of steps already convected for 'WKE'. Initialization sets StepCnt = 1
    INTEGER, INTENT(IN)                                            :: WkeSize       ! Number of total steps predefined for 'WKE'.
    DOUBLE PRECISION, INTENT (IN)                                              :: TStep        ! Time step for convection.
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)                                 :: V_inf

    !Output Variables
    type(Wake), INTENT(INOUT)     :: WKE

    !Local Variables
    INTEGER    :: nNodes, NodesStride
    INTEGER    :: nSgmnts, SgmntsStride
    INTEGER    :: Posic ! Identifies the Block of nodes/segments where the new set is to be inserted
    INTEGER    :: SgmntsInit  ! Sets the number of the initialized segments in the wake, the ones that coincide with the convection line.
    INTEGER    :: NodesInit   ! Sets the number of the initialized nodes in the wake, the ones that coincide with the convection line.
    INTEGER    :: i,j,k   ! Loop indexes


    ! Insertion and actualization parameters
    nNodes = WKE%nNodalPoints
    NodesStride = size(NewNodes)
    nSgmnts = WKE%nVortexSegments
    SgmntsStride = size(NewSegments)
    SgmntsInit = (SgmntsStride - 1)/2
    NodesInit = NodesStride   ! This variable is not extrictly necessary but is set to match the algorithm for Nodes and Sgmnts.

    if (StepCnt < WkeSize+1) then  !Covers the cases when there's still room in the wake for convection, no elimination needed.

        Wke%LLStride = Wke%LLStride - 1 !Identifies the last loaded stride position when there's still room in the wake (negative)
        !print *,'LLStride ',Wke%LLStride
        Posic = StepCnt

        DO i=0,Posic-1
            DO j=1,NodesStride ! Actualizes existent nodes positions
                WKE%NodalPointSet(i*NodesStride+j)%xyz = WKE%NodalPointSet(i*NodesStride+j)%xyz + (WKE%NodalPointSet(i*NodesStride+j)%InducedSpeed + V_inf) * TStep
            END DO
        END DO

        DO j=1,NodesStride !Inserts new block of NODES into the wake array which is not saturated (there are still room for adding nodes)
            WKE%NodalPointSet(Posic*NodesStride+j)%Label = NewNodes(j)%Label + NodesStride * Posic
            WKE%NodalPointSet(Posic*NodesStride+j)%xyz = NewNodes(j)%xyz
            WKE%NodalPointSet(Posic*NodesStride+j)%InducedSpeed = NewNodes(j)%InducedSpeed
        END DO

        DO k=1,SgmntsStride !Inserts new block of SEGMENTS into the wake array which is not saturated (there are still room for adding segments)
            WKE%VortexSegmentSet((Posic-1)*SgmntsStride+SgmntsInit+k)%Label = NewSegments(k)%Label + SgmntsStride * (Posic-1) + SgmntsInit
            WKE%VortexSegmentSet((Posic-1)*SgmntsStride+SgmntsInit+k)%Grid = NewSegments(k)%Grid
            WKE%VortexSegmentSet((Posic-1)*SgmntsStride+SgmntsInit+k)%SegmentNodes(1) = NewSegments(k)%SegmentNodes(1) + NodesStride * Posic !Makes the same increment as the one in nodes labels to reference the node ID according to the wake numeration
            WKE%VortexSegmentSet((Posic-1)*SgmntsStride+SgmntsInit+k)%SegmentNodes(2) = NewSegments(k)%SegmentNodes(2) + NodesStride * Posic !Makes the same increment as the one in nodes labels to reference the node ID according to the wake numeration
            WKE%VortexSegmentSet((Posic-1)*SgmntsStride+SgmntsInit+k)%SgmntCirc = NewSegments(k)%SgmntCirc
        END DO

        DO k=1,SgmntsInit  !Adds the circulation of the parallel segments of NewSegments to those for the previous loaded stride.
            WKE%VortexSegmentSet((Posic-1)*SgmntsStride+k)%SgmntCirc = WKE%VortexSegmentSet((Posic-1)*SgmntsStride+k)%SgmntCirc - NewSegments(SgmntsStride-SgmntsInit+k)%SgmntCirc !lllllllllllllllllllllllllllllllllllllllllllllllllllllll
        END DO

    else    ! Covers the cases when the wake is full and in order to add steps some blocks must be eliminated: (StepCnt>WkeSize)

        !Posic = mod(StepCnt,WkeSize + 1)
        Posic = mod(StepCnt-1,WkeSize)
        Wke%LLStride = Posic
        ! print *,'LLStrideAcc ',Wke%LLStride
        !print *, "POSIC: ", Posic
        DO i=0,Posic-1
            DO j=1,NodesStride ! Actualizes existing nodes positions and Labels for elements previous to the block insertion,except InitNodes
                WKE%NodalPointSet(i*NodesStride+NodesInit+j)%Label = WKE%NodalPointSet(i*NodesStride+NodesInit+j)%Label - NodesStride
                WKE%NodalPointSet(i*NodesStride+NodesInit+j)%xyz = WKE%NodalPointSet(i*NodesStride+NodesInit+j)%xyz + (WKE%NodalPointSet(i*NodesStride+NodesInit+j)%InducedSpeed + V_inf) * TStep
            END DO
        END DO

        !if (Posic<WkeSize) then
        DO i=0,Posic-1
            DO k=1,SgmntsStride ! Actualizes existing SEGMENTS Labels and Connectivities for elements previous to the block insertion, except InitSegments.
                WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%Label = WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%Label - SgmntsStride
                WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(1) = WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(1) - NodesStride
                WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(2) = WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(2) - NodesStride
            END DO
        END DO
        !end if

        !i = Posic
        !if (Posic<WkeSize) then
        DO k=1,NodesInit !Stores the last line of nodes into the 'init' positions for closing the last block of panels.
            WKE%NodalPointSet(k)%Label = WKE%NodalPointSet((Posic+1)*NodesStride+k)%Label  ! (Posic+1)= ...+NodesInit+...
            WKE%NodalPointSet(k)%xyz = WKE%NodalPointSet((Posic+1)*NodesStride+k)%xyz
            WKE%NodalPointSet(k)%InducedSpeed = WKE%NodalPointSet((Posic+1)*NodesStride+k)%InducedSpeed
        END DO

        DO k=1,NodesInit !Actualizes InitNodes lables and positions
            WKE%NodalPointSet(k)%Label = WKE%NodalPointSet(k)%Label - NodesStride
            WKE%NodalPointSet(k)%xyz = WKE%NodalPointSet(k)%xyz + (WKE%NodalPointSet(k)%InducedSpeed + V_inf) * TStep
        END DO

        DO j=1,NodesStride !Inserts new block of NODES into the wake array, which is saturated. Existing nodes are overwritten.
            WKE%NodalPointSet(Posic*NodesStride+NodesInit+j)%Label = nNodes - NodesStride + j
            WKE%NodalPointSet(Posic*NodesStride+NodesInit+j)%xyz = NewNodes(j)%xyz
            WKE%NodalPointSet(Posic*NodesStride+NodesInit+j)%InducedSpeed = NewNodes(j)%InducedSpeed
        END DO


        !if (Posic<WkeSize) then
        DO k=1,SgmntsInit !Stores the last line of parallel segments into the 'init' positions for closing the last block of panels.
            WKE%VortexSegmentSet(k)%Label = WKE%VortexSegmentSet((Posic+1)*SgmntsStride+k)%Label
            WKE%VortexSegmentSet(k)%Grid = WKE%VortexSegmentSet((Posic+1)*SgmntsStride+k)%Grid
            WKE%VortexSegmentSet(k)%SegmentNodes(1) = WKE%VortexSegmentSet((Posic+1)*SgmntsStride+k)%SegmentNodes(1)
            WKE%VortexSegmentSet(k)%SegmentNodes(2) = WKE%VortexSegmentSet((Posic+1)*SgmntsStride+k)%SegmentNodes(2)
            WKE%VortexSegmentSet(k)%SgmntCirc = WKE%VortexSegmentSet(k)%SgmntCirc - WKE%VortexSegmentSet((Posic+1)*SgmntsStride+k)%SgmntCirc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111111
        END DO

        DO k=1,SgmntsInit ! Actualizes the InitSegments Labels and Connectivities.
            WKE%VortexSegmentSet(k)%Label = WKE%VortexSegmentSet(k)%Label - SgmntsStride
            WKE%VortexSegmentSet(k)%SegmentNodes(1) = WKE%VortexSegmentSet(k)%SegmentNodes(1) - NodesStride
            WKE%VortexSegmentSet(k)%SegmentNodes(2) = WKE%VortexSegmentSet(k)%SegmentNodes(2) - NodesStride
        END DO

        DO k=1,SgmntsStride !Inserts new block of SEGMENTS into the wake array which is saturated. Existing nodes are overwritten
            WKE%VortexSegmentSet((Posic)*SgmntsStride+SgmntsInit+k)%Label = nSgmnts - SgmntsStride + k
            WKE%VortexSegmentSet((Posic)*SgmntsStride+SgmntsInit+k)%Grid = NewSegments(k)%Grid
            WKE%VortexSegmentSet((Posic)*SgmntsStride+SgmntsInit+k)%SegmentNodes(1) = NewSegments(k)%SegmentNodes(1) + NodesStride * WkeSize !Makes the same increment as the one in nodes labels to reference the node ID according to the wake numeration
            WKE%VortexSegmentSet((Posic)*SgmntsStride+SgmntsInit+k)%SegmentNodes(2) = NewSegments(k)%SegmentNodes(2) + NodesStride * WkeSize !Makes the same increment as the one in nodes labels to reference the node ID according to the wake numeration
            WKE%VortexSegmentSet((Posic)*SgmntsStride+SgmntsInit+k)%SgmntCirc = NewSegments(k)%SgmntCirc
        END DO

        if (Posic == 0) then
            DO k=1,SgmntsInit  !Adds the circulation of the parallel segments of NewSegments to those of the previous loaded stride.
                WKE%VortexSegmentSet(nSgmnts-SgmntsInit+k)%SgmntCirc = WKE%VortexSegmentSet(nSgmnts-SgmntsInit+k)%SgmntCirc - NewSegments(SgmntsStride-SgmntsInit+k)%SgmntCirc !Parallel new segments have the same circ with negative sign.
            END DO
        else
            DO k=1,SgmntsInit  !Adds the circulation of the parallel segments of NewSegments to those of the previous loaded stride.
                WKE%VortexSegmentSet((Posic)*SgmntsStride+k)%SgmntCirc = WKE%VortexSegmentSet((Posic)*SgmntsStride+k)%SgmntCirc - NewSegments(SgmntsStride-SgmntsInit+k)%SgmntCirc
            END DO
        end if
        !end if

        !DO k=1,SgmntsInit  !Adds the circulation of the parallel segments of NewSegments to those of the previous loaded stride.
        !    WKE%VortexSegmentSet((Posic)*SgmntsStride+k)%SgmntCirc = - NewSegments(SgmntsStride-SgmntsInit+k)%SgmntCirc
        !END DO

        DO i=Posic+1,WkeSize-1
            DO j=1,NodesStride ! Actualizes existing nodes positions and Labels for elements following to the block insertion.
                WKE%NodalPointSet(i*NodesStride+NodesInit+j)%Label = WKE%NodalPointSet(i*NodesStride+NodesInit+j)%Label - NodesStride
                WKE%NodalPointSet(i*NodesStride+NodesInit+j)%xyz = WKE%NodalPointSet(i*NodesStride+NodesInit+j)%xyz + (WKE%NodalPointSet(i*NodesStride+NodesInit+j)%InducedSpeed + V_inf) * TStep
            END DO
        END DO

        !if (Posic<WkeSize) then
        DO i=Posic+1,WkeSize-1
            DO k=1,SgmntsStride ! Actualizes existing SEGMENTS Labels and connectivities for elements following to the block insertion.
                WKE%VortexSegmentSet((i)*SgmntsStride+SgmntsInit+k)%Label = WKE%VortexSegmentSet((i)*SgmntsStride+SgmntsInit+k)%Label - SgmntsStride
                WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(1) = WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(1) - NodesStride
                WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(2) = WKE%VortexSegmentSet(i*SgmntsStride+SgmntsInit+k)%SegmentNodes(2) - NodesStride
            END DO
        END DO
        !end if

    end if

    ! Deallocation of auxilliary arrays
    deallocate(NewNodes)
    deallocate(NewSegments)

    END SUBROUTINE Wake_Actualization

    SUBROUTINE ConvectNodes_Adjust (WKE, GRD, ConvectionLine, StpCnt)
    ! Actualices the position of the last loaded nodes of wake according to the grid movement.
    ! The wake nodes of the convection line must move with the grid to avoid gaps.
    ! Wake-Grid-ConvectionLine relation is assumed.

    !Input Variables
    type(Convection_Line), INTENT(IN)  :: ConvectionLine
    type(Grid), INTENT(IN)             :: GRD
    INTEGER, INTENT(IN)                :: StpCnt

    !Output Variables
    type(Wake), INTENT(INOUT)       :: WKE

    !Local Variables
    INTEGER              :: NodeID    !Node label for identification
    INTEGER              :: Posic     !Identifies the Block of nodes/segments where the new set is to be inserted
    DOUBLE PRECISION, DIMENSION(3)   :: NewCoords !Variable extraction. Nodes coordinates
    INTEGER              :: NodesStride, NodesInit, WkeSize
    INTEGER              :: i         !Loop index

    NodesStride = WKE%NodesStride
    NodesInit = WKE%NodesStride
    WkeSize = WKE%WakeSize

    if (StpCnt < WkeSize+2) then  !Covers the cases when there's still room in the wake for convection, no elimination needed.

        Posic = StpCnt - 1 !Aims for the stride previously loaded

        DO i=1,ConvectionLine%nNodalPoints
            NodeID = ConvectionLine%NodalPointSet(i)
            !print*,'Nid = ', NodeID
            NewCoords =  GRD%NodalPointSet(NodeID)%xyz
            WKE%NodalPointSet(Posic*NodesStride+i)%xyz = NewCoords
            !print*,'Wid = ', WKE%NodalPointSet(Posic*NodesStride+i)%Label
        END DO

    else    ! Covers the cases when the wake is full and in order to add steps some blocks must be eliminated: (StepCnt>WkeSize)
        Posic = Wke%LLStride
        !print*,'LLstride = ', Wke%LLStride
        DO i=1,ConvectionLine%nNodalPoints
            NodeID = ConvectionLine%NodalPointSet(i)
            !print*,'Nid = ', NodeID
            NewCoords =  GRD%NodalPointSet(NodeID)%xyz
            WKE%NodalPointSet(Posic*NodesStride+NodesInit+i)%xyz = NewCoords
            !print*,'Wid = ', WKE%NodalPointSet(Posic*NodesStride+i)%Label
        END DO


    end if

    END SUBROUTINE ConvectNodes_Adjust


    SUBROUTINE Wake_Reading (WakeName, WakeLabel, WakeSize, Owner, ConvectLineID, ConvectLine_File, ExtraData)
    ! Read the wake's variables from the Body file. Couples with Wake allocation, initalization subroutines.
    ! Continues reading the "Body File" opened  in 'BodyRead' subroutine. No further input variables are needed.

    !Input Variables
    ! None

    !Output Variables
    CHARACTER(len=20),INTENT(OUT)        :: WakeName
    INTEGER,INTENT(OUT)                  :: WakeLabel
    INTEGER,INTENT(OUT)                  :: WakeSize
    INTEGER,INTENT(OUT)                  :: Owner            !Grid from which the wake is shed.
    INTEGER,INTENT(OUT)                  :: ConvectLineID
    CHARACTER(len=20),INTENT(OUT)        :: ConvectLine_File !Assumes CL is defined by nodes labels.
    CHARACTER(len=15),INTENT(OUT)        :: ExtraData

    !Local Variables
    INTEGER  :: ios  ! File handling check variable

    !Variables reading
    read (2,*,IOSTAT=ios) WakeLabel              !1
    print *, 'Wke Label..:', WakeLabel
    read (2,*,IOSTAT=ios) WakeName               !2
    print *, 'Wke Name..:', WakeName
    read (2,*,IOSTAT=ios) Owner                  !3
    print *, 'Wke Owner..:', Owner
    read (2,*,IOSTAT=ios) ConvectLineID          !4
    print *, 'Wke CL ID..:', ConvectLineID
    read (2,*,IOSTAT=ios) ConvectLine_File       !5
    print *, 'CL FILE..:', ConvectLine_File
    read (2,*,IOSTAT=ios) WakeSize               !6
    print *, 'WakeSize..:', WakeSize
    read (2,*,IOSTAT=ios) ExtraData

    pause

    END SUBROUTINE Wake_Reading


    SUBROUTINE Wake_Loading (WakeName, WakeLabel, WakeSize, ConvectLineID, ConvectLine_File, GridOwnr, WKE, ExtraData)
    !Loads a Wake into the Body structure using allocation and initialization subroutines.
    !Assumes that the Owner (Grid) is already loaded.
    !Includes the Convection Line processing.

    !Inputs Variables
    CHARACTER(len=20),INTENT(IN)        :: WakeName
    INTEGER,INTENT(IN)                  :: WakeLabel
    INTEGER,INTENT(IN)                  :: WakeSize
    INTEGER,INTENT(IN)                  :: ConvectLineID
    CHARACTER(len=20),INTENT(IN)        :: ConvectLine_File !Assumes CL is defined by nodes labels.
    CHARACTER(len=15),INTENT(IN)        :: ExtraData

    !Outputs Variables
    type(WAKE), INTENT(INOUT)           :: WKE
    type(Grid), INTENT(INOUT)           :: GridOwnr         !Grid from which the wake is shed. (Segment Condition application)

    !Local Variables
    INTEGER, ALLOCATABLE, DIMENSION(:) :: CL_VtxSgmnts ! Convection Line Segments Lables for CL loading.
    INTEGER                            :: nVtxSgmnts   ! Number of Convection Line Segments
    type(Convection_Line)              :: ConvectLine  ! Extraction Variable.

    !Convection Line processing
    print *, 'ConvectLine_File : ',ConvectLine_File
    call ConvecLine_Reading_N (ConvectLine_File, GridOwnr, CL_VtxSgmnts, nVtxSgmnts)
    call ConvecLine_Loading (GridOwnr, CL_VtxSgmnts, nVtxSgmnts, ConvectLineID, ConvectLine)

    !Grid Segment Condition --> Classifies the grid segments according to convective characteristics
    print *, 'Sgmnt_condition : ', GridOwnr%GridLabel, ConvectLine%CLineLabel
    call Segment_Condition (GridOwnr, ConvectLine)

    !Wake Loading
    WKE%WakeName = WakeName
    WKE%WakeLabel = WakeLabel
    WKE%WakeSize = WakeSize
    WKE%ExtraData = ExtraData

    print *, 'Wake_Alloc : '
    call Wake_Allocation (WKE, ConvectLine, WakeSize)
    print *,'Wake_Init : '
    call Wake_Initialization (WKE, ConvectLine, GridOwnr)
    print *, 'Wake_Load_END : '
    END SUBROUTINE Wake_Loading


    END MODULE classWake