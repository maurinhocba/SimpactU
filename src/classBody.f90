    MODULE classBody

    ! BODY abstract data type definition

    !Represents the widest class of the library.

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza


    USE classGrid!, only: Grid
    USE classBC_Grid!, only: BC_Grid
    USE classWake!, only: Wake
    USE classVortex_Segment, only: Vortex_Segment
    USE classNodal_Point, only: Nodal_Point
    USE classPanel, only: Panel
    USE classExtra_Point!, only: Extra_Point

    IMPLICIT NONE
    PUBLIC


    TYPE, Public :: Body

        CHARACTER(len=20)                               :: BodyName
        INTEGER                                         :: BodyLabel ! Body ID, sets the reference for its components
        INTEGER                                         :: nGrids   ! Number of Grids of the Body
        TYPE(Grid), ALLOCATABLE, DIMENSION(:)           :: Grids ! Grid Set of the Body
        INTEGER                                         :: nBCGrids !Number of BCGrids of the body
        INTEGER, ALLOCATABLE, DIMENSION(:)              :: GlblGridLabels
        INTEGER, ALLOCATABLE, DIMENSION(:)              :: GlblBCGridLabels
        INTEGER, ALLOCATABLE, DIMENSION(:)              :: GlblWakeLabels
        TYPE(BC_Grid), ALLOCATABLE, DIMENSION (:)       :: BCGrids ! BCGrid Set of the Body
        INTEGER                                         :: nWakes   !Number of Wakes of the Body
        TYPE(Wake), ALLOCATABLE, DIMENSION(:)           :: Wakes ! Wake Set of the Body
        INTEGER                                         :: nExtraPoints ! Number of ExtraPoints of the Body
        TYPE(Extra_Point), ALLOCATABLE, DIMENSION(:)    :: ExtraPoints
        CHARACTER(len=20)                               :: extraData  ! for misc purposes

    CONTAINS ! methods' names
    PROCEDURE, Public, NoPass :: Body_ReadLoad
    PROCEDURE, Public, NoPass :: BodySet_AllocGen

    END TYPE Body

    CONTAINS ! ==================================================================

    SUBROUTINE Body_ReadLoad (Body_FileName, V_inf, BDY)
    !Reads BodyFile and sets the Bodies Variables into BodySet.
    !Needs linking with Grid, BCGrid, and Wake loading/constructing subroutines.
    !Continues on reading the file opened by 'MainVars_Read' subroutine. No in_file positioning is needed.

    !Input variables
    CHARACTER(len=20), INTENT(IN)  :: Body_FileName ! Main variables file name
    !INTEGER, INTENT(IN)            :: NumberOfBodies, NumberOfGrids, NumberOfBCGrids  ! Number of structures
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: V_inf

    !Output variables
    type(BODY), INTENT(INOUT) :: BDY

    !Local variables
    INTEGER  :: ios  ! File handling check variable
    INTEGER  :: i, j, k, l, m    ! Loop indexes
    CHARACTER(len=20) :: ExtrPntFile
    CHARACTER(len=15) :: ExtraData


    CHARACTER (len=15)                    :: GridName
    INTEGER                               :: GridLabel
    INTEGER                               :: NumberOfNodes, NumberOfPanels!, NumberOfSegments
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: NodesArray
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: IenArray
    INTEGER                               :: NumberOfWakes
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: WakesLabels
    INTEGER                               :: KinCondition

    CHARACTER(len=15)                     :: BCGridName
    INTEGER                               :: BCGridLabel
    INTEGER                               :: BCNumberOfNodes
    INTEGER                               :: BCNumberOfPanels
    DOUBLE PRECISION                                  :: TranspirationVel
    INTEGER                               :: TrnspCondition
    CHARACTER(len=20)                     :: TranspVel_FileName !>>> Not processed yet !!!!!!!!!!!!!!!!!!!!!11
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: BCNodesArray
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: BCIenArray

    CHARACTER(len=20)        :: WakeName
    INTEGER                  :: WakeLabel
    INTEGER                  :: WakeSize
    INTEGER                  :: ConvectLineID
    INTEGER                  :: Owner
    INTEGER                  :: OwAux    ! Aux variable for grid position in array extraction
    CHARACTER(len=20)        :: ConvectLine_File
    !type(Grid)               :: GridOwnr

    !BodyFile open
    open (2, FILE=Body_FileName, IOSTAT=ios, ACTION='read')

    ! Opening status check
    if (ios /= 0) then

        print * , 'File opening error - Body File: ',Body_FileName
        pause
    else
        print * ,'Body File initialized... '
    end if

    read (2,*,IOSTAT=ios) BDY%BodyName   !1
    read (2,*,IOSTAT=ios) BDY%BodyLabel  !2
    read (2,*,IOSTAT=ios) BDY%nGrids     !3

    allocate (BDY%GlblGridLabels(BDY%nGrids)) ! Array of Grids Global Variables allocation
    allocate (BDY%Grids(BDY%nGrids))

    read (2,*,IOSTAT=ios) BDY%nBCGrids   !4

    allocate (BDY%GlblBCGridLabels(BDY%nBCGrids)) ! Array of BCGrids Global Variables allocation
    allocate (BDY%BCGrids(BDY%nBCGrids))

    read (2,*,IOSTAT=ios) BDY%nWakes   !5
    print *, 'nWakes-> : ',BDY%nWakes, '--'


    allocate (BDY%GlblWakeLabels(BDY%nWakes)) ! Array of Wakes Global Variables allocation
    allocate (BDY%Wakes(BDY%nWakes))

    read (2,*,IOSTAT=ios) BDY%extraData   !6

    read (2,*,IOSTAT=ios) BDY%GlblGridLabels  !7
    print *, 'GlblGridLabels : ',BDY%GlblGridLabels, '--'

    read (2,*,IOSTAT=ios) BDY%GlblBCGridLabels  !8
    print *, 'GlblBCGridLabels : ',BDY%GlblBCGridLabels, '--'

    read (2,*,IOSTAT=ios) BDY%GlblWakeLabels  !9
    print *, 'GlblWakeLabels : ',BDY%GlblWakeLabels, '--'

    read (2,*,IOSTAT=ios) BDY%nExtraPoints !10

    allocate (BDY%ExtraPoints(BDY%nExtraPoints)) ! Array of ExtraPoints allocation

    read (2,*,IOSTAT=ios) ExtrPntFile !11   - Extra Point Coordinates file name. Must be set as argument of Extra Point loading subrutine for the body.

    !call EXTRA POINT LOADING SUBROUTINE: ExtraPointReading, ExtraPointLoading
    print *, 'BDY%nExtPoints ',BDY%nExtraPoints
    if (BDY%nExtraPoints > 0) then
        call ExtraPoint_ReadLoad (ExtrPntFile, BDY%nExtraPoints, BDY%ExtraPoints)
    end if


    !Check print
    print *, 'BodyVars:  ', BDY%BodyName,BDY%BodyLabel,BDY%nGrids,BDY%nBCGrids,BDY%nWakes,BDY%GlblGridLabels(:),BDY%GlblBCGridLabels(:),BDY%GlblWakeLabels(:),BDY%nExtraPoints,ExtrPntFile(:)

    DO j=1, BDY%nGrids
        call GridReading (GridName, GridLabel, NumberOfNodes, NumberOfPanels, NumberOfWakes, WakesLabels, KinCondition, NodesArray, IENArray, ExtraData)
        call GridLoading (BDY%Grids(j), GridName, GridLabel, NumberOfNodes, NumberOfPanels, NumberOfWakes, WakesLabels, KinCondition, NodesArray, IenArray, V_inf, ExtraData)
    END DO

    DO k=1, BDY%nBCGrids
        print *, 'BC_Grids'
        call BC_GridReading (BCGridName, BCGridLabel, BCNumberOfNodes, BCNumberOfPanels, TrnspCondition, KinCondition, TranspirationVel, TranspVel_FileName, BCNodesArray, BCIENArray, ExtraData)
        call BC_GridLoading (BDY%BCGrids(k), BCGridName, BCGridLabel, BCNumberOfNodes, BCNumberOfPanels, KinCondition, TrnspCondition, TranspirationVel, TranspVel_FileName, BCNodesArray, BCIENArray, V_inf, ExtraData)
        print *, 'BCGrid k : ',k
    END DO

    print *, 'BDY%nWakes ',BDY%nWakes
    DO l=1, BDY%nWakes
        print *, 'Wakes'
        call Wake_Reading (WakeName, WakeLabel, WakeSize, Owner, ConvectLineID, ConvectLine_File, ExtraData)
        !Owner Search
        DO i=1, BDY%nGrids
            if (BDY%Grids(i)%GridLabel == Owner) then
                print *, 'Owner.....GrdLabel i',BDY%Grids(i)%GridLabel,i
                OwAux = i
                exit
            end if
        END DO
        call Wake_Loading (WakeName, WakeLabel, WakeSize, ConvectLineID, ConvectLine_File, BDY%Grids(OwAux), BDY%Wakes(l), ExtraData)
    END DO



    pause
    !File Close
    close (UNIT=2)

    END SUBROUTINE Body_ReadLoad


    SUBROUTINE BodySet_AllocGen (NumberOfBodies, V_inf, BodyFile_list, Body_Set)   ! Set allocation and Generation
    !Allocates and generates the body set for the simulation.

    !Input Variables
    INTEGER, INTENT(IN)                        :: NumberOfBodies ! Number of bodies of the simulation.
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)             :: V_inf          ! Free stream velocity  .
    CHARACTER(len=20), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: BodyFile_list ! List of body files for data input.

    !Output Variables
    type(BODY),ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: Body_Set ! Array set variable for allocation

    !Local Variables
    INTEGER :: i ! Loop Index
    print *, 'BodyFile_list..: ',BodyFile_list
    pause

    ALLOCATE (Body_Set (NumberOfBodies))

    DO i=1,NumberOfBodies
        call Body_ReadLoad (BodyFile_list(i), V_inf, Body_Set(i))
    END DO

    END SUBROUTINE BodySet_AllocGen

    END MODULE classBody

