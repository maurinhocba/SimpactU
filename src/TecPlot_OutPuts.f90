MODULE TecPlot_OutPuts

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    ! Contains subroutines for TecPlot format output files generation.

    USE classGrid
    USE classVortex_Segment
    USE classNodal_Point
    USE classPanel
    USE classWake
    USE classConvection_Line
    USE classBody

    IMPLICIT NONE
    PUBLIC


    CONTAINS

SUBROUTINE TPOP_FileInit (FileUnit, OP_FileName)
    ! Generates a TecPlot format output generic file given by 'FileUnit'. Sets the headings for further printing.
    ! TPOP: TecPlot OutPut

    !Input Variables
    INTEGER, INTENT(IN)            :: FileUnit !Destination File
    CHARACTER(len=20), INTENT (IN) :: OP_FileName    ! Debug file name
    !Local Variables
    INTEGER            :: ios  ! File opening check
    CHARACTER (len=15) :: CaseName  ! Name for job identification

    Open (UNIT=FileUnit, FILE=OP_FileName, Status='New', Action='Write', Iostat=ios)  ! ------------FILE NAME MUST BE PROCESSED!!!!!
    ! Opening status check
    if (ios /= 0) then
        print * , 'File opening error - OutPut File: ',OP_FileName
        pause
    else
        print * ,'TecPlot Output File initialized... '
    end if


    CaseName = 'TRIAL CASE 1'

    ! Header writting
    write (FileUnit,'(A20,A,A7,I10,A1)') 'TITLE = "CASE: ', CaseName,'"'
    write (FileUnit,'(A34,///)') 'VARIABLES = "X" "Y" "Z" "DeltaCP"'! "n(X)" "n(Y)" "n(Z)" "Cp"

    Close (UNIT = FileUnit)

    END SUBROUTINE TPOP_FileInit

    SUBROUTINE TPOP_Step (BodySet,nBodies,StpCnt,FileUnit)
    !Gathers the Tecplot OutPut subroutines for each step. Includes Girds, BcGrids and Wakes.

    !Input Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(IN) :: BodySet
    INTEGER, INTENT(IN)     :: FileUnit !Destination File
    INTEGER, INTENT(IN)     :: nBodies  ! Number of Bodies of the Set
    INTEGER, INTENT(IN)     :: StpCnt ! Step Counter = Solution Time
    !Output Variables
    !NONE

    !Local Variables
    INTEGER :: i,j,k,l !Loop Index
    INTEGER :: StrandID, SIDoffset

    SIDoffset = 0

    DO i=1,nBodies !  First loop on Bodies
        DO j=1,BodySet(i)%nGrids !  First loop on Grids and BCGrids
            StrandID=SIDoffset + j
            call Grid_TPOP (BodySet(i)%Grids(j),StpCnt,StrandID,FileUnit)
        END DO
        DO k=1,BodySet(i)%nBCGrids
            StrandID=SIDoffset + k + BodySet(i)%nGrids
            call Grid_TPOP (BodySet(i)%BCGrids(k),StpCnt,StrandID,FileUnit)
        END DO
        DO l=1,BodySet(i)%nWakes
            StrandID=SIDoffset + l + BodySet(i)%nGrids + BodySet(i)%nBCGrids
            call Wake_TPOP (BodySet(i)%Wakes(l),StrandID,StpCnt,BodySet(i)%Wakes(l)%WakeSize,BodySet(i)%Wakes(l)%ConvecLine,FileUnit)
        END DO
        !SIDoffset = BodySet(i)%nGrids+BodySet(i)%nBCGrids+BodySet(i)%nWakes
        SIDoffset = StrandID
    END DO


 END SUBROUTINE TPOP_Step

SUBROUTINE ExtPts_TPOP_FileInit (FileUnit, OP_FileName)
    ! Generates a TecPlot format output generic file FOR EXTRA POINTS processing given by 'FileUnit'. Sets the headings for further printing.
    ! TPOP: TecPlot OutPut

    !Input Variables
    INTEGER, INTENT(IN)            :: FileUnit !Destination File
    CHARACTER(len=20), INTENT (IN) :: OP_FileName    ! File name
    !Local Variables
    INTEGER            :: ios  ! File opening check
    CHARACTER (len=20) :: CaseName  ! Name for job identification

    Open (UNIT=FileUnit, FILE=OP_FileName, Status='New', Action='Write', Iostat=ios)  ! ------------FILE NAME MUST BE PROCESSED!!!!!
    ! Opening status check
    if (ios /= 0) then
        print * , 'File opening error - OutPut File: ',OP_FileName
        pause
    else
        print * ,'TecPlot ExtraPoints Output File initialized... '
    end if


    CaseName = 'Extra_Points Output'

    ! Header writting
    write (FileUnit,'(A20,A,A7,I10,A1)') 'TITLE = "CASE: ', CaseName,'"'
    write (FileUnit,'(A42,//)') 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "|v|"'

    Close (UNIT = FileUnit)

END SUBROUTINE ExtPts_TPOP_FileInit

SUBROUTINE ExtPnt_TPOP (ExtPoint_Set, StrandID, StpCnt, FileUnit)
    ! Generates a TecPlot format output file for Extra Points.
    ! TPOP: TecPlot OutPut

    !Input Variables
    INTEGER, INTENT(IN)            :: FileUnit !Destination File
    !CHARACTER(len=20), INTENT (IN) :: OP_FileName    ! File name
    type(Extra_Point), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ExtPoint_Set
    INTEGER, INTENT (IN)									 :: StrandID, StpCnt

    !Local Variables
    INTEGER            :: ios  ! File opening check
    CHARACTER (len=30) :: CaseName  ! Name for job identification
    !INTEGER            :: FileUnit ! Destination File
    INTEGER            :: NumberOfExtPnt
    INTEGER			   :: i
    CHARACTER(len=20)					 :: ExtPnt_OPfile    ! file name

    NumberOfExtPnt = size(ExtPoint_Set)
    ExtPnt_OPfile = 'ExtPnt_OPF.dat'

    Open (UNIT=FileUnit, FILE=ExtPnt_OPfile, Status='Old', Position='Append', Action='Write', Iostat=ios)  ! ------------FILE NAME MUST BE PROCESSED!!!!! ,Status='New'
    ! Opening status check
    if (ios /= 0) then
        print * , 'File opening error - OutPut File: ',ExtPnt_OPfile
        pause
    else
        print * ,'TecPlot EXTRA POINT Output File initialized... '
    end if

    !Zone definition
    write (FileUnit,'(A26,I7,A2)') 'ZONE T= "ExtraPoints  ID: ',StrandID,'",'
    write (FileUnit,'(A16,I7)') 'STRANDID= 201'!, StrandID
    write (FileUnit,'(A16,I7)') 'SOLUTIONTIME= ', StpCnt
    write (FileUnit,'(A20)') 'DATAPACKING = BLOCK'
    write (FileUnit,'(A4,I15)') 'I = ', NumberOfExtPnt
    write (FileUnit,'(A48)') 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)  '


    ! Extra Point Coordinates
    DO i=1, NumberOfExtPnt
        write (FileUnit,'(F16.4,$)') ExtPoint_Set(i)%xyz(1)
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if
    END DO
    write (FileUnit,'(/)') !Line space

    DO i=1, NumberOfExtPnt
        write (FileUnit,'(F16.4,$)') ExtPoint_Set(i)%xyz(2)
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if
    END DO
    write (FileUnit,'(/)') !Line space

    DO i=1, NumberOfExtPnt
        write (FileUnit,'(F16.4,$)') ExtPoint_Set(i)%xyz(3)
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if
    END DO
    write (FileUnit,'(/)') !Line space

    !ExtrPoints Speed
    DO i=1, NumberOfExtPnt
        write (FileUnit,'(F15.5,A2,$)') ExtPoint_Set(i)%Velocity(1),'  '
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if

    END DO

    DO i=1, NumberOfExtPnt
        write (FileUnit,'(F15.5,A2,$)') ExtPoint_Set(i)%Velocity(2), '  '
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if

    END DO

    DO i=1, NumberOfExtPnt
        write (FileUnit,'(F15.5,A2,$)')  ExtPoint_Set(i)%Velocity(3), '  '
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if

    END DO

    DO i=1, NumberOfExtPnt
        write (FileUnit,'(F15.5,A2,$)')  ExtPoint_Set(i)%VelocityMag, '  '
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if

    END DO

    END SUBROUTINE ExtPnt_TPOP

    SUBROUTINE Grid_TPOP (GRD,StpCnt,StrandID,FileUnit)

    ! Generates a TecPlot format output file for a Grid or a BC_GRID.
    ! USES POLYMORPHIC VARIABLES SO AS TO PROCESS GRIDS AND BC_GRIDS.
    ! TPOP: TecPlot OutPut

    !Input Variables
    CLASS(Grid),INTENT (IN) :: GRD !POLYMORPHIC VARIABLE
    !type(Grid), INTENT (IN) :: GRD
    INTEGER, INTENT(IN)     :: FileUnit !Destination File
    INTEGER, INTENT(IN)     :: StrandID ! ID for multiple zone output.
    INTEGER, INTENT(IN)     :: StpCnt ! Step Counter = Solution Time

    !Local Variables
    INTEGER  :: i,j  !loop index
    INTEGER  :: ios !File opening check

    !File Opening
    Open (UNIT=FileUnit, FILE='A_OPTP.dat', Status='Old', Position='Append', Action='Write', Iostat=ios)
    print *, 'ios: ', ios  ! File opening check
    print *, 'GRID NAME: ', GRD%GridName

    !Zone definition
    write (FileUnit,'(A8,A,A5,I7,A2)') 'ZONE T="GRID:',GRD%GridName,' ID: ',StrandID,'",'
    write (FileUnit,'(A12,I7)') 'STRANDID= ', StrandID
    write (FileUnit,'(A16,I7)') 'SOLUTIONTIME= ', StpCnt

    if (GRD%PanelSet(1)%Panel_Type == 3) then
        write (FileUnit,'(A27)') 'ZONETYPE = FETriangle'
    else
        write (FileUnit,'(A27)') 'ZONETYPE = FEQuadrilateral'
    end if

    write (FileUnit,'(A4,I5,A4,I5,A1)') 'E = ',GRD%nPanels,',N = ',GRD%nNodalPoints,','
    write (FileUnit,'(A31)') 'VARLOCATION=([4]=CELLCENTERED)'
    write (FileUnit,'(A10,A15)') 'C = ',GRD%ExtraData
    write (FileUnit,'(A20)') 'DATAPACKING = BLOCK'
    write (FileUnit,'(A33)') 'DT=(SINGLE SINGLE SINGLE DOUBLE)'


    ! Nodes Coordinates
    DO i=1, GRD%nNodalPoints
        write (FileUnit,'(F16.4,$)') GRD%NodalPointSet(i)%xyz(1)
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if
    END DO
    write (FileUnit,'(/)') !Line space

    DO i=1, GRD%nNodalPoints
        write (FileUnit,'(F16.4,$)') GRD%NodalPointSet(i)%xyz(2)
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if
    END DO
    write (FileUnit,'(/)') !Line space

    DO i=1, GRD%nNodalPoints
        write (FileUnit,'(F16.4,$)') GRD%NodalPointSet(i)%xyz(3)
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if
    END DO
    write (FileUnit,'(/)') !Line space

    !DeltaCP
    DO i=1, GRD%nPanels
        write (FileUnit,'(F25.5,A2,$)') GRD%PanelSet(i)%CtrlPoint%Delta_Cp, '  '
        if (mod(i,50) == 0) then
            write (FileUnit,'(/)') !Line space, for line length reduction
        end if

    END DO

    write (FileUnit,'(/)') !Line space

    !PanelNodes Conectivities     - Supports only 4-nodes quadrilaterals (PanleType=1)
    !if (StpCnt == 1) then
    DO i=1, GRD%nPanels

        DO j=1,GRD%PanelSet(i)%nPanelNodes
            write (FileUnit,'(I10,$)') GRD%PanelSet(i)%PanelNodes(j)

            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO

    END DO
    write (FileUnit,'(/)') !Line space
    !end if

END SUBROUTINE Grid_TPOP

SUBROUTINE Wake_TPOP (WKE,StrandID,StpCnt,WkeSize,ConvecLine,FileUnit)

    ! Generates a TecPlot format output file for a Wake.
    ! TPOP: TecPlot OutPut

    !Input Variables
    type(Wake), INTENT (IN)          :: WKE
    type(Convection_Line),INTENT(IN) :: ConvecLine
    INTEGER, INTENT(IN)              :: FileUnit !Destination File
    INTEGER, INTENT(IN)              :: StrandID ! ID for multiple zone output.
    INTEGER, INTENT(IN)              :: StpCnt ! Step Counter = Solution Time
    INTEGER, INTENT(IN)              :: WkeSize

    !Local Variables
    INTEGER  :: i  !loop index
    INTEGER  :: ios !File opening check
    INTEGER  :: NodesStride ! Number of convection line nodes. Used for first node location for printing.
    INTEGER  :: SgmntsStride ! Number of segments added per step. Used for first node location for printing.
    INTEGER  :: SgmntsInit   ! Number of convection line segments. Used for first segment location for printing.
    INTEGER  :: NodesInit   ! Number of convection line nodes. Used for first node location for printing.
    INTEGER  :: LoopInit ! First node coordinate position in the array.

    NodesStride = ConvecLine%nNodalPoints
    SgmntsStride = 2*ConvecLine%nVortexSegments+1
    SgmntsInit = ConvecLine%nVortexSegments
    NodesInit = NodesStride

    !File Opening
    Open (UNIT=FileUnit, FILE='A_OPTP.dat', Status='Old', Position='Append', Action='Write', Iostat=ios)
    print *, 'ios: ', ios  ! File opening check


    !Zone definition
    write (FileUnit,'(A8,A,A5,I7,A2)') 'ZONE T="WAKE:',WKE%WakeName,' ID: ',StrandID,'",'
    write (FileUnit,'(A10,I7)') 'STRANDID= ', StrandID
    write (FileUnit,'(A14,I7)') 'SOLUTIONTIME= ', StpCnt
    write (FileUnit,'(A32)') 'ZONETYPE = FEQuadrilateral,'
    if (StpCnt < WkeSize) then
        write (FileUnit,'(A4,I5,A4,I5,A1)') 'E = ',StpCnt*SgmntsStride+SgmntsInit,',N = ',(StpCnt+1)*NodesStride,','
    else
        write (FileUnit,'(A4,I5,A4,I5,A1)') 'E = ',WKE%nVortexSegments,',N = ',WKE%nNodalPoints,','
    end if
    write (FileUnit,'(A19)') 'PASSIVEVARLIST=[4]'
    write (FileUnit,'(A20)') 'DATAPACKING = BLOCK'
    write (FileUnit,'(A27)') 'DT=(SINGLE SINGLE SINGLE)'
    write (FileUnit,'(A10,A15)') 'C = ', WKE%ExtraData

    ! Nodes Coordinates
    if (StpCnt <= WkeSize) then
        DO i=1, StpCnt*NodesStride + NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(1), '  '

            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if
        END DO
        write (FileUnit,'(/)') !Line space

        DO i=1, StpCnt*NodesStride + NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(2), '  '

            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        DO i=1, StpCnt*NodesStride + NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(3), '  '

            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if
        END DO
        write (FileUnit,'(/)') !Line space

    else
        DO i=1,NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(1), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        LoopInit = mod(StpCnt-1,WkeSize) + 1 ! Position of last stride loaded ("Posic") plus one stride.
        DO i=LoopInit*NodesStride+NodesInit+1,WKE%nNodalPoints
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(1), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        DO i=NodesInit+1, LoopInit*NodesStride+NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(1), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        DO i=1,NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(2), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        LoopInit = mod(StpCnt-1,WkeSize) + 1 ! Position of last stride loaded ("Posic") plus one stride.
        DO i=LoopInit*NodesStride+NodesInit+1,WKE%nNodalPoints
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(2), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        DO i=NodesInit+1, LoopInit*NodesStride+NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(2), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        DO i=1,NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(3), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        LoopInit = mod(StpCnt-1,WkeSize) + 1 ! Position of last stride loaded ("Posic") plus one stride.
        DO i=LoopInit*NodesStride+NodesInit+1,WKE%nNodalPoints
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(3), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space

        DO i=NodesInit+1, LoopInit*NodesStride+NodesInit
            write (FileUnit,'(F20.5,A2,$)') WKE%NodalPointSet(i)%xyz(3), '  '
            if (mod(i,50) == 0) then
                write (FileUnit,'(/)') !Line space, for line length reduction
            end if

        END DO
        write (FileUnit,'(/)') !Line space
    end if

    write (FileUnit,'(/)') !Line space

    !PanelNodes Conectivities     - Supports only 4-nodes quadrilaterals (PanleType=1)
    if (StpCnt<= WkeSize) then
        DO i=1, StpCnt*SgmntsStride+SgmntsInit
            write (FileUnit,'(4I12)') WKE%VortexSegmentSet(i)%SegmentNodes(1),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(1)
        END DO
    else

        DO i=1,SgmntsInit
            write (FileUnit,'(4I12)') WKE%VortexSegmentSet(i)%SegmentNodes(1),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(1)
        END DO

        LoopInit = mod(StpCnt-1,WkeSize) + 1 ! Position of last stride loaded ("Posic") plus one stride.
        DO i=LoopInit*SgmntsStride+SgmntsInit+1,WKE%nVortexSegments
            write (FileUnit,'(4I12)') WKE%VortexSegmentSet(i)%SegmentNodes(1),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(1)
        END DO
        DO i=SgmntsInit+1, LoopInit*SgmntsStride+SgmntsInit
            write (FileUnit,'(4I12)') WKE%VortexSegmentSet(i)%SegmentNodes(1),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(2),WKE%VortexSegmentSet(i)%SegmentNodes(1)
        END DO

    end if

END SUBROUTINE Wake_TPOP

END MODULE TecPlot_OutPuts


