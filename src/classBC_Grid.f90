    MODULE classBC_Grid

    ! "BC_Grid" extended abstract data type definition
    ! Specialization of "Grid"

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE classGrid!, only: grid

    IMPLICIT NONE
    PUBLIC


    TYPE, EXTENDS(Grid), Public :: BC_Grid  !Extended type of "Grid"
        ! Attributes:

        INTEGER ::  TrnspCondition   ! Transpiration velocity flag: No penetration (0); Constant (1); Variable (2). For variable trnsp vel ControlPoint TrnspVel must be activated.
        DOUBLE PRECISION    ::  TranspirationVel !Transpiration Velocity (normal to surface). =0 if TrnspCondition=0; =constant if TrnspCondition=1 ; not used if TrnspCondition=2 (read value from file and load into CP)
        CHARACTER(len=20) :: TranspVel_FileName !>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!11

    CONTAINS ! methods' names
    PROCEDURE, Public, NoPass :: BC_GridLoading
    !PROCEDURE, Public, NoPass :: BC_GridPrinting
    PROCEDURE, Public, NoPass :: BC_GridReading
    PROCEDURE, Public, NoPass :: BC_CPTrnspVelLoading

    END TYPE BC_Grid

    CONTAINS ! ==================================================================

    SUBROUTINE BC_GridReading (BCGridName, BCGridLabel, NumberOfNodes, NumberOfPanels, TrnspCondition, KinCondition, TranspirationVel, TranspVel_FileName, NodesArray, IENArray, ExtraData)
    !Reads the body file for extracting BCGrid variables. Couples with 'BCGridLoading' subroutine.
    !Continues reading the "Body File" opened  in 'BodyRead' subroutine. No further input variables are needed.

    !Input Variables
    !None

    !Output Variables
    CHARACTER(len=15), INTENT(OUT)                     :: BCGridName
    INTEGER, INTENT(OUT)                               :: BCGridLabel
    INTEGER, INTENT(OUT)                               :: NumberOfNodes
    INTEGER, INTENT(OUT)                               :: NumberOfPanels
    DOUBLE PRECISION, INTENT(OUT)                                  :: TranspirationVel
    INTEGER, INTENT(OUT)                               :: TrnspCondition
    INTEGER, INTENT(OUT)                               :: KinCondition
    CHARACTER(len=20), INTENT(OUT)                     :: TranspVel_FileName !>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!11
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT)    :: NodesArray
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT) :: IenArray
    CHARACTER(len=15),INTENT(OUT)                      :: ExtraData

    !Local Variables
    INTEGER  :: ios  ! File opening check variable
    CHARACTER(len=20) :: NA_FileName ! Nodes coordinates file name.
    CHARACTER(len=20) :: IEN_FileName! Connectivities array file name.
    INTEGER :: PnlType !>>>> Not used yet
    INTEGER :: NCPvelFlag !>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!!11
    CHARACTER(len=20) :: NCPvel_FileName!>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!!11

    !File Reading - BC_Grid Variables

    read (2,*,IOSTAT=ios) BCGridLabel       !1
    read (2,*,IOSTAT=ios) BCGridName        !2
    read (2,*,IOSTAT=ios) NumberOfNodes     !3
    read (2,*,IOSTAT=ios) NumberOfPanels    !4

    read (2,*,IOSTAT=ios) NA_FileName       !7
    read (2,*,IOSTAT=ios) PnlType           !8   >>>> Not used yet
    read (2,*,IOSTAT=ios) IEN_FileName      !9
    read (2,*,IOSTAT=ios) KinCondition      !10
    read (2,*,IOSTAT=ios) NCPvel_FileName   !11 !!!!!!!!!!!!!!!!!!!! Process
    read (2,*,IOSTAT=ios) TrnspCondition    !12
    if (TrnspCondition==1) then
        read (2,*,IOSTAT=ios) TranspirationVel  !13 a
        print *,'TrnspVel val: ',TranspirationVel
    else
        read (2,*,IOSTAT=ios) TranspVel_FileName !13 b
        print *,'TrnspVel file: ',TranspVel_FileName
    end if
    read (2,*,IOSTAT=ios) ExtraData

    !Check print
    print *, 'BCGridVars',BCGridName, BCGridLabel, NumberOfNodes, NumberOfPanels, TrnspCondition, IEN_FileName, NA_FileName


    ! CALL NODES_ARRAY AND IEN_ARRAY READING
    call IEN_Array_R(IEN_FileName, NumberOfPanels, IenArray)
    call Nodes_Array_R(NA_FileName,NumberOfNodes,NodesArray)

    END SUBROUTINE BC_GridReading


    SUBROUTINE BC_GridLoading (BCGRD, BCGridName, BCGridLabel, NumberOfNodes, NumberOfPanels, KinCondition, TrnspCondition, TranspirationVel, TranspVel_FileName, NodesArray, IENArray, V_inf, ExtraData)
    !Constructs a BC_GRID class from the 'BC_GridGenerator/Reading' output data. Uses as inputs the coordinates, connectivities,
    !and lists arrays of the structure. Also uses P_S_N module for further data constructs.

    !Input Variables
    CHARACTER(len=15), INTENT(IN)                     :: BCGridName
    INTEGER, INTENT(IN)                               :: BCGridLabel
    INTEGER, INTENT(IN)                               :: NumberOfNodes
    INTEGER, INTENT(IN)                               :: NumberOfPanels
    DOUBLE PRECISION, INTENT(IN)                                  :: TranspirationVel
    INTEGER, INTENT(IN)                               :: TrnspCondition
    INTEGER, INTENT(IN)                               :: KinCondition
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)                    :: V_inf
    CHARACTER(len=20), INTENT(IN)                     :: TranspVel_FileName !>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!11
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:), INTENT(IN)    :: NodesArray
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(IN) :: IenArray
    CHARACTER(len=15), INTENT(IN)                     :: ExtraData

    !Output Variables
    type(BC_Grid), INTENT (INOUT) :: BCGRD

    !Local Variables
    INTEGER :: i !loop index
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: IESArray
    INTEGER                               :: NumberOfSegments
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ISNArray
    INTEGER, ALLOCATABLE, DIMENSION (:)   :: ShPanels_Links, ShPanels_Storage
    INTEGER, DIMENSION (:), ALLOCATABLE   :: PanelType


    call Panel_Type (IENArray, PanelType)
    call Segments (BCGridLabel, IENArray, PanelType, NumberOfSegments, ISNArray, IESArray,ShPanels_Links, ShPanels_Storage)

    print *, 'BCName -------'
    !Grid identification
    BCGRD%GridName = BCGridName
    BCGRD%GridLabel = BCGridLabel
    print *, 'BCConfig -------'
    !Grid configuration
    BCGRD%KinCondition = KinCondition
    BCGRD%TrnspCondition = TrnspCondition
    print *, 'BC_NPs -------'
    !NodalPoints Loading
    BCGRD%nNodalPoints = NumberOfNodes
    call NPSet_Allocation (NumberOfNodes, BCGRD%NodalPointSet)
    call NPSet_Generation (BCGRD%NodalPointSet, NodesArray)
    print *, 'BCVSs -------'
    !VortexSegments Loading
    BCGRD%nVortexSegments= NumberOfSegments
    call VSSet_Allocation (NumberOfSegments, BCGRD%VortexSegmentSet)
    call VSSet_Generation (BCGridLabel, IENArray, BCGRD%VortexSegmentSet, NumberOfSegments, ISNArray, IESArray, ShPanels_Links, ShPanels_Storage, PanelType)
    print *, 'BCPnls -------'
    !Panels Loading
    BCGRD%nPanels= NumberOfPanels
    call PanelSet_Allocation (NumberOfPanels, BCGRD%PanelSet)
    call PanelSet_Generation (BCGridLabel, IENArray, IESArray, NodesArray, PanelType, BCGRD%PanelSet, V_inf)
    print *, 'BCCPs -------'
    !CALL CP TRANSP VELOCITY LOADING (FILE/CONSTANT/NO PENETRATION)----------------------------------<<<<<<
    call BC_CPTrnspVelLoading (BCGRD,TranspVel_FileName,TranspirationVel,TrnspCondition)
    print *, 'BCsgmnt coef -------'
    !Segments Orientation
    do i=1, BCGRD%nPanels
        call Sgmnt_Coef (BCGRD%PanelSet(i), BCGRD%VortexSegmentSet, BCGRD%NodalPointSet)
    end do
    print *, 'BCsgmnt cond -------'
    !Segments Condition: Replaces Grids' Segment_Condition subroutine.
    DO i=1,BCGRD%nVortexSegments
        BCGRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 1.0
        BCGRD%VortexSegmentSet(i)%Sgmnt_Cond(2) = 0.0

        !Edge Identification
        if (BCGRD%VortexSegmentSet(i)%nSharingPanels == 1) then
            BCGRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 2.0
        end if
    END DO
    print *, 'BCpnl coef -------'
    ! Panels Orientations
    call Panel_Coef(BCGRD)

    BCGRD%ExtraData = ExtraData

    print *, 'Name: ', BCGRD%GridName
    print *, 'Label: ', BCGRD%GridLabel
    print *, 'Is rigid?: ', BCGRD%KinCondition
    print *, 'Number Of Nodes: ', BCGRD%nNodalPoints
    print *, 'Number Of Segments: ', BCGRD%nVortexSegments
    print *, 'Number Of Panels: ', BCGRD%nPanels
    !pause

    END SUBROUTINE BC_GridLoading

    !SUBROUTINE BC_Panel_Coef (BCGRD)
    !    ! Generates an array of coeficient (1,-1) for the panels sharing the segment ('SharingPanels') indicating the sign of the contribution
    !    ! of the ring vorticity of the panel to the circulation of the segment. The structure is applied on the panel list of each segment ('SharingPanels')
    !    ! adding a negative sign to the panel label if it's contribution to the segment circulation is negative. For example, the list SharingPanels=[253,-24,16,-89] indicates that the panel
    !    ! #24 and #89 substract their vorticity to the segment circulation while the panels #253 and #16 add to it.
    !    ! The procedure works on a grid since it modifies its segment set using data from its panel set.
    !
    !!Input Variable
    !type(BC_Grid),INTENT (INOUT)    :: BCGRD
    !
    !!Local Variables
    !INTEGER     :: i,j,k   !Loop indexes
    !INTEGER     :: nShPnls, nPnlSgmnts  ! Number of sharing panels, number of panel segments. Extraction variables
    !INTEGER     :: PanelID ! Panel Label. Extraction variable
    !INTEGER     :: SGN     ! Sign of the label. Extraction variable
    !
    !DO i=1,BCGRD%nVortexSegments
    !    nShPnls = BCGRD%VortexSegmentSet(i)%nSharingPanels
    !
    !    DO j=1,nShPnls
    !        PanelID = BCGRD%VortexSegmentSet(i)%SharingPanels(j)
    !        nPnlSgmnts = BCGRD%PanelSet(PanelID)%nPanelSegments
    !        DO k=1,nPnlSgmnts
    !            if (abs(BCGRD%PanelSet(PanelID)%PanelSegments(k))==BCGRD%VortexSegmentSet(i)%Label) then
    !                SGN = sign(1,BCGRD%PanelSet(PanelID)%PanelSegments(k))
    !            end if
    !        END DO
    !        BCGRD%VortexSegmentSet(i)%SharingPanels(j) = SGN * BCGRD%VortexSegmentSet(i)%SharingPanels(j)
    !    END DO
    !END DO
    !
    !END SUBROUTINE BC_Panel_Coef

    SUBROUTINE BC_CPTrnspVelLoading (BCGRD,TranspVel_FileName,TranspirationVel,TrnspCondition)
    !Reads Transpiration velocity file and loads the values on each BC_Grid CP (CP%TrnspVelocity)
    !Could be reused if the Transpiration Velocity is time dependant, appending data without closing the file.

    !Input Variables
    DOUBLE PRECISION, INTENT(IN)                      :: TranspirationVel ! Always normal to surface. Is considered the proyected value.
    INTEGER, INTENT(IN)                   :: TrnspCondition   ! Transp Velocity condition. Defines the read/load process to apply.
    CHARACTER(len=20), INTENT(IN)         :: TranspVel_FileName ! File containig the values of each CP's transpiration velocity in correspondance to the panel order.

    !Output Variables
    type(BC_Grid), INTENT(INOUT) :: BCGRD

    !Local Variables
    INTEGER :: i !loop index
    INTEGER :: ios !File handling check variable

    if (TrnspCondition == 2) then
        open (7, FILE=TranspVel_FileName, IOSTAT=ios, ACTION='read')
        if (ios /= 0) then
            print * , 'File opening error - TranspVel File: ',TranspVel_FileName
            pause
        end if

        DO i=1,BCGRD%nPanels
            read (7,*,IOSTAT=ios) BCGRD%PanelSet(i)%CtrlPoint%TrnspVelocity
        END DO

        close(UNIT=7) ! Avoid closing if the Transpiration velocity is time dependant and more data will be read from the file.

    else if (TrnspCondition == 1) then
        DO i=1,BCGRD%nPanels
            BCGRD%PanelSet(i)%CtrlPoint%TrnspVelocity = TranspirationVel
        END DO
    else
        DO i=1,BCGRD%nPanels
            BCGRD%PanelSet(i)%CtrlPoint%TrnspVelocity = 0.0
        END DO
    end if

    END SUBROUTINE BC_CPTrnspVelLoading

    END MODULE classBC_Grid
