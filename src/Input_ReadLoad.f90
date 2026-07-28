MODULE Input_ReadLoad
    
    !Gathers subroutines for Reading Input Files and Loading Variables into some structures.
    
  !USE classVortex_Segment
  !USE classNodal_Point
  !USE classPanel
  !!USE classGrid
  !USE classWake
  !USE classConvection_Line
  !!USE classBody
  
  IMPLICIT NONE
  PUBLIC
  
    CONTAINS

SUBROUTINE MainVars_Read (MVars_FileName, NSteps, TStep, VinfFile, OPfile, RSfile, DBfile, VinfFlag, OPflag, DBflag, RSflag, RS_time, V_inf, NumberOfBodies, WakeSize, NumberOfGrids, NumberOfBCGrids, NumberOfWakes, BodyFile_list)
!Reads MainInput File and sets the "Main Variables" for the program excecution.

!Input Variables
CHARACTER(len=20), INTENT(IN) :: MVars_FileName

!Output Variables
INTEGER, INTENT(OUT)                  :: NSteps  ! Number of Steps of the simulation
DOUBLE PRECISION, INTENT(OUT)                     :: TStep   ! Time duration of the step, time step.
INTEGER, INTENT(OUT)                  :: NumberOfBodies, NumberOfGrids, NumberOfBCGrids, NumberOfWakes ! Number of structures of the simulation: bodies, grids, bc_grids, wakes.
INTEGER, INTENT(OUT)                  :: WakeSize ! General number of steps for wakes. Different sizes for each wake can be set in each wake file.
DOUBLE PRECISION, DIMENSION(3), INTENT(OUT)       :: V_inf ! Free stream speed. Initial value or constant value. If it's time dependant must be read from corresponding file.
CHARACTER(len=20), INTENT(OUT)        :: VinfFile, OPfile, RSfile, DBfile ! Data files names variables.
INTEGER, INTENT(OUT)                  :: VinfFlag, OPflag, DBflag, RSflag ! Conditional reading flags for variables (Vinf, OutPut, DeBug, ReStart)
DOUBLE PRECISION, INTENT(OUT)                     :: RS_time ! Restart time, sets the CPU time for restart saving.
!type(BODY), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: BodySet ! Set of bodies of the simulation.
CHARACTER(len=20), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: BodyFile_list ! List of body files for data input.

!Local Variables
INTEGER  :: ios  ! File opening check variable
INTEGER  :: i    ! Loop index

!Open File
open (1, FILE=MVars_FileName, IOSTAT=ios, ACTION='read')

! Opening status check
    if (ios /= 0) then
    
        print * , 'File opening error - MainVars File: ',MVars_FileName
        pause
    
    end if
    
!File Reading - Main Variables
        
    read (1,*,IOSTAT=ios) NSteps          !1
    read (1,*,IOSTAT=ios) TStep           !2
    read (1,*,IOSTAT=ios) VinfFlag        !3
    read (1,*,IOSTAT=ios) VinfFile        !4
    read (1,*,IOSTAT=ios) V_inf            !5 x3
    read (1,*,IOSTAT=ios) OPflag          !6
    read (1,*,IOSTAT=ios) OPfile          !7
    read (1,*,IOSTAT=ios) DBflag          !8
    read (1,*,IOSTAT=ios) DBfile          !9
    read (1,*,IOSTAT=ios) RSflag          !10
    read (1,*,IOSTAT=ios) RS_time         !11
    read (1,*,IOSTAT=ios) RSfile          !12
    read (1,*,IOSTAT=ios) NumberOfBodies  !13
    
    allocate(BodyFile_list(NumberOfBodies))  ! List of bodies file names allocation
    
    read (1,*,IOSTAT=ios) WakeSize        !14
    read (1,*,IOSTAT=ios) NumberOfGrids   !15
    read (1,*,IOSTAT=ios) NumberOfBCGrids !16
    read (1,*,IOSTAT=ios) NumberOfWakes   !17
    read (1,*,IOSTAT=ios) BodyFile_list   !18


!BodySet allocation
!allocate (BodySet (NumberOfBodies)) 

 !Check print
 print *, NSteps,TStep,VinfFlag,VinfFile,V_inf,OPflag,OPfile,DBflag,DBfile,RSflag,RS_time,RSfile,NumberOfBodies,WakeSize,NumberOfGrids,NumberOfBCGrids,NumberOfWakes,BodyFile_list(:)
 !pause
 
 !do i=1,numberofbodies
 !    call body_read (bodyfile_list(i), v_inf, bodyset(i))
 !end do
 

!File Close
 close (UNIT=1)      
 
END SUBROUTINE MainVars_Read
           
SUBROUTINE Nodes_Array_R(NA_FileName,NumberOfNodes,NodesArray)
!Reads a Coordinates file for Nodes of a grid and constucts a NodesArray.

!Input Variables
CHARACTER(len=20), INTENT(IN)   :: NA_FileName ! Nodes coordinates file name.
INTEGER,INTENT(IN)      :: NumberOfNodes      ! Number of nodes in the file.

!Output Variables
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: NodesArray ! Array of coordinates "Nodes Array"

!Local Variables
INTEGER    :: ios !File check variable
INTEGER    :: i !Loop index
DOUBLE PRECISION, DIMENSION(3) :: NodeInLine ! File line extraction variable
print * , 'NodesArray File: ',NA_FileName
!Open File
open (1, FILE=NA_FileName, IOSTAT=ios, ACTION='read')

! Opening status check
    if (ios /= 0) then
    
        print * , 'File opening error - NodesArray File: ',NA_FileName
        pause
    
    end if

!Allocate NodesArray
    allocate (NodesArray (3, NumberOfNodes))

!File reading - Coordinates
 do i=1,NumberOfNodes
        
        read (1,*,IOSTAT=ios) NodeInLine
        NodesArray(1,i) = NodeInLine(1)
        NodesArray(2,i) = NodeInLine(2)
        NodesArray(3,i) = NodeInLine(3)
        
 end do
 
 !Check print
 do i=1,NumberOfNodes
        print *,NodesArray(:,i),'/'
 end do
 print * , 'NodesArray File: ',NA_FileName
 !File Close
 close (UNIT=1)
    
END SUBROUTINE Nodes_Array_R

SUBROUTINE IEN_Array_R(IEN_FileName, NumberOfPanels, IEN_Array)
    !Reads a Connectivities file for Panels of a grid and constucts a IENArray.

!Input Variables
CHARACTER(len=20), INTENT(IN)   :: IEN_FileName ! Panels Connectivities file name.
INTEGER,INTENT(IN)      :: NumberOfPanels      ! Number of Panels in the file.

!Output Variables
INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT) :: IEN_Array

!Local Variables
INTEGER    :: ios !File check variable
INTEGER    :: i,j !Loop index
INTEGER, DIMENSION(8) :: IENinLine ! File line extraction variable

!Open File
open (1, FILE=IEN_FileName, IOSTAT=ios, ACTION='read')

! Opening status check
    if (ios /= 0) then
    
        print * , 'File opening error - IENArray File: ',IEN_FileName
        pause
    
    end if

!Allocate IENArray
    allocate (IEN_Array (8, NumberOfPanels))
    !pause
!File reading - Connectivities
    IENinLine = 0
 do i=1,NumberOfPanels
        
        read (1,*,IOSTAT=ios) IENinLine
        do j=1,8

             IEN_Array(j,i) = IENinLine(j)
        end do
                
 end do
 
 !Check print
 do i=1,NumberOfPanels
        print *,IEN_Array(:,i),'/'
 end do
 !pause
 !File Close
 close (UNIT=1)

END SUBROUTINE IEN_Array_R
    
    END MODULE INPUT_ReadLoad