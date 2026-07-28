    MODULE classGrid

    ! GRID abstract data type definition

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE classVortex_Segment
    USE classNodal_Point
    USE classPanel
    USE P_S_N
    USE Input_ReadLoad

    IMPLICIT NONE
    PUBLIC


    TYPE, Public :: Grid

        CHARACTER (len=15)                              :: GridName
        INTEGER                                         :: GridLabel ! Grid ID. Sets the reference for components (nodes, panels, vtxsgments...)
        CHARACTER (len=15)                              :: ExtraData  ! for miscellaneous purposes
        INTEGER                                         :: KinCondition ! Kinematic condition of the grid: free(0), rigid(1), fixed(2).
        INTEGER                                         :: nWakes ! Number of wakes shedding from the grid.
        INTEGER, ALLOCATABLE, DIMENSION(:)              :: WakesLabels ! Labels of the nWakes wakes shedding from the grid.
        INTEGER                                         :: nNodalPoints ! number of nodal points
        type(NODAL_POINT), ALLOCATABLE, DIMENSION(:)    :: NodalPointSet ! (Set) Nodal points of the grid
        INTEGER                                         :: nVortexSegments ! number of vortex segments
        type(VORTEX_SEGMENT), ALLOCATABLE, DIMENSION(:) :: VortexSegmentSet ! (Set) Vortex segments of the grid
        INTEGER                                         :: nPanels ! number of panels
        type(PANEL), ALLOCATABLE, DIMENSION(:)          :: PanelSet  ! (Set) Panels of the grid
        INTEGER                                         :: AMatLoc ! Global AeroMatrix position index.

    CONTAINS ! methods' names

    PROCEDURE, Public, NoPass :: GridLoading
    PROCEDURE, Public, NoPass :: GridPrinting
    PROCEDURE, Public, NoPass :: InGrid_Xa
    PROCEDURE, Public, NoPass :: Panel_Coef
    PROCEDURE, Public, NoPass :: RingVrtct_2_SgmntCirc
    PROCEDURE, Public, NoPass :: PanelSet_RingVtct_Loading
    PROCEDURE, Public, NoPass :: GridReading

    END TYPE Grid

    CONTAINS ! ==================================================================

    SUBROUTINE GridReading (GridName, GridLabel, NumberOfNodes, NumberOfPanels, NumberOfWakes, WakesLabels, KinCondition, NodesArray, IENArray, ExtraData)
    !Reads the body file for extracting grid variables. Couples with 'GridLoading' subroutine.
    !Continues reading the "Body File" opened  in 'BodyRead' subroutine. No further input variables are needed.

    !Input Variables
    !None

    !Output Variables
    CHARACTER(len=15), INTENT(OUT)                     :: GridName
    INTEGER, INTENT(OUT)                               :: GridLabel
    INTEGER, INTENT(OUT)                               :: NumberOfNodes
    INTEGER, INTENT(OUT)                               :: NumberOfPanels
    INTEGER, INTENT(OUT)                               :: NumberOfWakes
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: WakesLabels
    INTEGER, INTENT(OUT)                               :: KinCondition
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT)    :: NodesArray
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT) :: IenArray
    CHARACTER(len=15), INTENT(OUT)                     :: ExtraData

    !Local Variables
    INTEGER  :: ios  ! File opening check variable
    CHARACTER(len=20) :: NA_FileName ! Nodes coordinates file name.
    CHARACTER(len=20) :: IEN_FileName! Connectivities array file name.
    INTEGER :: PnlType !>>>> Not used yet
    INTEGER :: NCPvelFlag !>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!!11
    CHARACTER(len=20) :: NCPvel_FileName!>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!!11

    !File Reading - Main Variables

    read (2,*,IOSTAT=ios) GridLabel             !1
    read (2,*,IOSTAT=ios) GridName          !2
    read (2,*,IOSTAT=ios) NumberOfNodes     !3
    read (2,*,IOSTAT=ios) NumberOfPanels    !4
    read (2,*,IOSTAT=ios) NumberOfWakes     !5

    allocate (WakesLabels(NumberOfWakes))
    read (2,*,IOSTAT=ios) WakesLabels       !6

    read (2,*,IOSTAT=ios) NA_FileName       !7
    read (2,*,IOSTAT=ios) PnlType           !8   >>>> Not used yet
    read (2,*,IOSTAT=ios) KinCondition      !9
    read (2,*,IOSTAT=ios) NCPvelFlag        !10 !!!!!!!!!!!!!!!!!!!! Process
    read (2,*,IOSTAT=ios) IEN_FileName      !11
    read (2,*,IOSTAT=ios) NCPvel_FileName   !12 !!!!!!!!!!!!!!!!!!!! Process
    read (2,*,IOSTAT=ios) ExtraData
    !Check print
    print *, 'GridVars:  ',GridName, GridLabel, NumberOfNodes, NumberOfPanels, NumberOfWakes, WakesLabels, KinCondition, IEN_FileName, NA_FileName


    ! CALL NODES_ARRAY AND IEN_ARRAY READING
    call IEN_Array_R(IEN_FileName, NumberOfPanels, IenArray)
    call Nodes_Array_R(NA_FileName,NumberOfNodes,NodesArray)

    END SUBROUTINE GridReading

    SUBROUTINE GridLoading (GRD, GridName, GridLabel, NumberOfNodes, NumberOfPanels, NumberOfWakes, WakesLabels, KinCondition, NodesArray, IenArray, V_inf, ExtraData)
    !Constructs a GRID object from the 'GridGenerator/Reading' output data. Uses as inputs the coordinates, connectivities,
    !and lists arrays of the structure. Also uses P_S_N module for further data constructs.

    !Input Variables
    CHARACTER (len=15), INTENT(IN)                    :: GridName
    INTEGER                                           :: GridLabel
    INTEGER, INTENT(IN)                               :: NumberOfNodes, NumberOfPanels!, NumberOfSegments
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:), INTENT(IN)    :: NodesArray
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(IN) :: IenArray
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)                    :: V_inf
    INTEGER, INTENT(IN)                               :: NumberOfWakes
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)    :: WakesLabels
    INTEGER, INTENT(IN)                               :: KinCondition
    CHARACTER(len=15), INTENT(IN)                     :: ExtraData

    !Output Variables
    type(Grid), INTENT (INOUT) :: GRD

    !Local Variables
    INTEGER :: i !loop index
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: IESArray
    INTEGER                               :: NumberOfSegments
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ISNArray
    INTEGER, ALLOCATABLE, DIMENSION (:)   :: ShPanels_Links, ShPanels_Storage
    INTEGER, DIMENSION (:), ALLOCATABLE   :: PanelType


    call Panel_Type (IENArray, PanelType)
    call Segments (GridLabel, IENArray, PanelType, NumberOfSegments, ISNArray, IESArray,ShPanels_Links, ShPanels_Storage)


    !Grid identification
    GRD%GridName = GridName
    GRD%GridLabel = GridLabel

    !Grid configuration
    GRD%nWakes = NumberOfWakes
    GRD%WakesLabels = WakesLabels
    GRD%KinCondition = KinCondition
    GRD%ExtraData = ExtraData

    !NodalPoints Loading
    GRD%nNodalPoints = NumberOfNodes
    call NPSet_Allocation (NumberOfNodes, GRD%NodalPointSet)
    call NPSet_Generation (GRD%NodalPointSet, NodesArray)

    !VortexSegments Loading
    GRD%nVortexSegments= NumberOfSegments
    call VSSet_Allocation (NumberOfSegments, GRD%VortexSegmentSet)
    call VSSet_Generation (GridLabel, IENArray, GRD%VortexSegmentSet, NumberOfSegments, ISNArray, IESArray, ShPanels_Links, ShPanels_Storage, PanelType)

    !Panels Loading
    GRD%nPanels= NumberOfPanels
    call PanelSet_Allocation (NumberOfPanels, GRD%PanelSet)
    call PanelSet_Generation (GridLabel, IENArray, IESArray, NodesArray, PanelType, GRD%PanelSet, V_inf)

    !Segments Orientation
    do i=1, GRD%nPanels
        call Sgmnt_Coef (GRD%PanelSet(i), GRD%VortexSegmentSet, GRD%NodalPointSet)
    end do

    !Segment Length
    call VSSet_Length (GRD%VortexSegmentSet, GRD%NodalPointSet)

    ! Panels Orientations
    call Panel_Coef(GRD)


    print *, 'Name: ', GRD%GridName
    print *, 'Label: ', GRD%GridLabel
    print *, 'Is rigid?: ', GRD%KinCondition
    print *, 'Number Of Nodes: ', GRD%nNodalPoints
    print *, 'Number Of Segments: ', GRD%nVortexSegments
    print *, 'Number Of Panels: ', GRD%nPanels
    pause

    END SUBROUTINE GridLoading

    SUBROUTINE PanelSet_RingVtct_Loading (GRD, G)
    ! Stores the ring vorticity calculated by the SOLVER in te current step, into the grid panels.
    ! Also stores the previous step rig vorticity of the panels for DeltaCP calculation. 'RngVtctyPrevStep'
    ! It is required for further calculations over the grid in the aerodynamic steps.
    !USES POLYMORPHIC VARIABLE SO AS TO PROCESS GRIDS AND BC_GRIDS

    !Input Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: G

    !Output Variables
    CLASS(GRID), INTENT(INOUT)  :: GRD
    !type(GRID), INTENT(INOUT)  :: GRD

    !Local Variables
    INTEGER :: i !Loop index

    DO i=1, GRD%nPanels
        !print*,'GPrev = ', GRD%PanelSet(i)%RingVorticity
        GRD%PanelSet(i)%RngVtctyPrvStp = GRD%PanelSet(i)%RingVorticity  ! Moves previous step value.
        GRD%PanelSet(i)%RingVorticity = G(i)                            ! Stores current step value.
    END DO

    END SUBROUTINE PanelSet_RingVtct_Loading

    SUBROUTINE GridPrinting (GRD, FileUnit)
    !Prints a Grid into a file. Requires the file to be opened and available for writing.

    !Input Variables
    INTEGER, INTENT(IN)    :: FileUnit
    type(Grid), INTENT(IN) :: GRD

    !Local Variables


    !Grid General Data Printing

    write (FileUnit,'(A19,/)') 'GRID General Data: '
    write (FileUnit,10) 'Name: ', GRD%GridName
    write (FileUnit,11) 'Label: ', GRD%GridLabel
    !write (FileUnit,12) 'Extra Data: ', GRD%ExtraData  ----> Extra Data is not defined---
    write (FileUnit,13) 'Is rigid?: ', GRD%KinCondition
    write (FileUnit,14) 'Number Of Nodes: ', GRD%nNodalPoints
    write (FileUnit,15) 'Number Of Segments: ', GRD%nVortexSegments
    write (FileUnit,16) 'Number Of Panels: ', GRD%nPanels
    write (FileUnit,'(/)')

    ! Nodal Point Set Printing
    call NPSet_Printing (GRD%NodalPointSet, FileUnit)

    ! Vortex Segments Set Printing
    call VSSet_Printing (GRD%VortexSegmentSet, FileUnit)

    ! Panel Set Printing
    call PanelSet_Printing (GRD%PanelSet, FileUnit)

    !FORMATS
10  FORMAT (A6,A)
11  FORMAT (A7,I10)
    !12 FORMAT (A12,A15) ----> Extra Data is not defined---
13  FORMAT (A11,L1)
14  FORMAT (A17,I9)
15  FORMAT (A17,I9)
16  FORMAT (A17,I9)



    END SUBROUTINE GridPrinting

    SUBROUTINE InGrid_Xa (PanelNodes, NodalPointSet, PanelType ,Xa)
    ! Calculates a matrix with the coordinates of each node of a panel as columns. Uses strcutures of a Grid (sets).

    !Input Variables
    INTEGER, DIMENSION (:), INTENT(IN)                       :: PanelNodes ! Panel Nodes ID
    type(NODAL_POINT), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: NodalPointSet  ! Set of Nodes of the grid.
    INTEGER, INTENT (IN)                                     :: PanelType  ! Sets the configuration of the Xa matrix.

    !Output Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT) :: Xa     ! Output matrix

    !Local Variables
    INTEGER :: nNodes, nID ! Number of nodes of the panel - ID of each node
    INTEGER :: i ! Loop index

    nNodes = size(PanelNodes)

    ALLOCATE (Xa(3,8))
    Xa = 0.d0 ! Initialize with zeros.

    if (PanelType == 3) then
        do i=1,3    !CASE 3.
            Xa(:,i) =  NodalPointSet(PanelNodes(i))%xyz
        end do
    else
        do i=1,4   !CASE 1. Applies for all except for CASE 3.
            Xa(:,i) =  NodalPointSet(PanelNodes(i))%xyz
        end do
    end if

    select case (PanelType)

    case (51)
        Xa(:,5) =  NodalPointSet(PanelNodes(5))%xyz
    case (52)
        Xa(:,6) =  NodalPointSet(PanelNodes(6))%xyz
    case (53)
        Xa(:,7) =  NodalPointSet(PanelNodes(7))%xyz
    case (54)
        Xa(:,8) =  NodalPointSet(PanelNodes(8))%xyz
    case (61)
        Xa(:,5) =  NodalPointSet(PanelNodes(5))%xyz
        Xa(:,6) =  NodalPointSet(PanelNodes(6))%xyz
    case (62)
        Xa(:,6) =  NodalPointSet(PanelNodes(6))%xyz
        Xa(:,7) =  NodalPointSet(PanelNodes(7))%xyz
    case (63)
        Xa(:,7) =  NodalPointSet(PanelNodes(7))%xyz
        Xa(:,8) =  NodalPointSet(PanelNodes(8))%xyz
    case (64)
        Xa(:,5) =  NodalPointSet(PanelNodes(5))%xyz
        Xa(:,8) =  NodalPointSet(PanelNodes(8))%xyz
    case (71)
        Xa(:,5) =  NodalPointSet(PanelNodes(5))%xyz
        Xa(:,6) =  NodalPointSet(PanelNodes(6))%xyz
        Xa(:,7) =  NodalPointSet(PanelNodes(7))%xyz
    case (72)
        Xa(:,6) =  NodalPointSet(PanelNodes(6))%xyz
        Xa(:,7) =  NodalPointSet(PanelNodes(7))%xyz
        Xa(:,8) =  NodalPointSet(PanelNodes(8))%xyz
    case (73)
        Xa(:,5) =  NodalPointSet(PanelNodes(5))%xyz
        Xa(:,7) =  NodalPointSet(PanelNodes(7))%xyz
        Xa(:,8) =  NodalPointSet(PanelNodes(8))%xyz
    case (74)
        Xa(:,5) =  NodalPointSet(PanelNodes(5))%xyz
        Xa(:,6) =  NodalPointSet(PanelNodes(6))%xyz
        Xa(:,8) =  NodalPointSet(PanelNodes(8))%xyz
    case (8)
        do i=5,8
            Xa(:,i) =  NodalPointSet(PanelNodes(i))%xyz
        end do

    end select

    END SUBROUTINE InGrid_Xa

    SUBROUTINE Panel_Coef (GRD)
    ! Generates an array of coeficient (1,-1) for the panels sharing the segment ('SharingPanels') indicating the sign of the contribution
    ! of the ring vorticity of the panel to the circulation of the segment. The structure is applied on the panel list of each segment ('SharingPanels')
    ! adding a negative sign to the panel label if it's contribution to the segment circulation is negative. For example, the list SharingPanels=[253,-24,16,-89] indicates that the panel
    ! #24 and #89 substract their vorticity to the segment circulation while the panels #253 and #16 add to it.
    ! The procedure works on a grid since it modifies its segment set using data from its panel set.
    ! USES POLYMORPHIC VARIABLES SO AS TO PROCESS GRIDS AND BC_GRIDS.

    !Input Variable
    CLASS(Grid),INTENT (INOUT)    :: GRD

    !Local Variables
    INTEGER     :: i,j,k   !Loop indexes
    INTEGER     :: nShPnls, nPnlSgmnts  ! Number of sharing panels, number of panel segments. Extraction variables
    INTEGER     :: PanelID ! Panel Label. Extraction variable
    INTEGER     :: SGN     ! Sign of the label. Extraction variable

    !pause
    DO i=1,GRD%nVortexSegments
        nShPnls = GRD%VortexSegmentSet(i)%nSharingPanels
        !print *, 'BCnShpnls: ',GRD%VortexSegmentSet(i)%nSharingPanels

        DO j=1,nShPnls
            PanelID = GRD%VortexSegmentSet(i)%SharingPanels(j)
            nPnlSgmnts = GRD%PanelSet(PanelID)%nPanelSegments
            !print *, 'PanelID , nPnlSgmnts: ',PanelID,GRD%PanelSet(PanelID)%nPanelSegments
            DO k=1,nPnlSgmnts
                if (abs(GRD%PanelSet(PanelID)%PanelSegments(k))==GRD%VortexSegmentSet(i)%Label) then
                    SGN = sign(1,GRD%PanelSet(PanelID)%PanelSegments(k))
                end if
            END DO
            GRD%VortexSegmentSet(i)%SharingPanels(j) = SGN * GRD%VortexSegmentSet(i)%SharingPanels(j)
        END DO
    END DO
    !print *,'end pnlCoef ---'
    END SUBROUTINE Panel_Coef

    SUBROUTINE RingVrtct_2_SgmntCirc (GRD)
    ! Calculates Vrotex Segment circulation from surrounding panels ring vorticity for a grid.
    ! Asumes 'Panel_Coef' subroutine previous excecution.
    ! Operates in a Grid/BC_Grid, setting its segments circulation from its Panel Set data.
    !USES POLYMORPHIC VARIABLE SO AS TO PROCESS GRIDS AND BC_GRIDS

    !Input Variables
    CLASS(Grid), INTENT(INOUT)    :: GRD
    !type(Grid), INTENT(INOUT)    :: GRD

    !Local Variables
    INTEGER   :: i,j        ! Loop indexes
    INTEGER   :: PanelID    ! Label of the panel. Extraction variable
    INTEGER   :: nShPnls    ! Number of sharing panels. Extraction variable
    DOUBLE PRECISION      :: RingVtct   ! Ring vorticity of the panel. Extraction Variable
    DOUBLE PRECISION      ::  b         ! Sgmnt Condition coeficient extraction variable
    DOUBLE PRECISION      :: RingVtctPS ! Ring vorticity of the panel at previous time step. Extraction Variable
    DOUBLE PRECISION      :: Gamma      ! Circulation of the segment. Extraction Variable
    INTEGER   :: SGN        ! Sign of the label. Extraction variable

    DO i=1,GRD%nVortexSegments
        nShPnls = GRD%VortexSegmentSet(i)%nSharingPanels
        Gamma = 0.0
        !b = GRD%VortexSegmentSet(i)%Sgmnt_Cond(2) ! Segment condition coeficient "b"********************

        DO j=1,nShPnls
            PanelID = abs(GRD%VortexSegmentSet(i)%SharingPanels(j))
            RingVtct = GRD%PanelSet(PanelID)%RingVorticity
            !RingVtctPS = GRD%PanelSet(PanelID)%RngVtctyPrvStp!*****************************

            ! Gc = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SgmntCirc ! Actual Sgmnt Circulation
            ! Gp = sign(1,Pnl%PanelSegments(i))*Pnl%RngVtctyPrvStp    ! Previous Sgmnt Circulation

            SGN = sign(1,GRD%VortexSegmentSet(i)%SharingPanels(j))

            Gamma = Gamma + SGN * RingVtct
            !Gamma = Gamma + SGN * (RingVtct - b*RingVtctPS)

        END DO

        GRD%VortexSegmentSet(i)%SgmntCirc = Gamma

    END DO

    END SUBROUTINE RingVrtct_2_SgmntCirc


    END MODULE classGrid
