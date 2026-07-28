    program Main

    USE ifport
    USE classNodal_Point
    USE classVortex_Segment
    USE classPanel
    USE classGrid
    USE classWake
    USE classBody
    USE P_S_N
    USE AeroMatrix
    USE TecPlot_OutPuts
    USE classConvection_Line
    USE STEP
    USE Aero_Loads
    USE Input_ReadLoad

    IMPLICIT NONE

    !Main Variables for the simulation
    INTEGER                                      :: NumberOfBodies, NumberOfGrids, NumberOfBCGrids, NumberOfWakes ! Number of structures of the simulation: bodies, grids, bc_grids, wakes.
    INTEGER                                      :: WakeSize ! General number of steps for wakes. Different sizes for each wake can be set in each wake file.
    DOUBLE PRECISION, DIMENSION(3)               :: V_inf ! Free stream speed. Initial value or constant value. If it's time dependant must be read from corresponding file.
    CHARACTER(len=20)                            :: VinfFile, OPfile, RSfile, DBfile, ExtPnt_OPfile ! Data files names variables.
    INTEGER                                      :: VinfFlag, OPflag, DBflag, RSflag ! Conditional reading flags for variables (Vinf, OutPut, DeBug, ReStart)
    DOUBLE PRECISION                             :: RS_time ! Restart time, sets the CPU time for restart saving.
    type(BODY), ALLOCATABLE, DIMENSION(:)        :: Body_Set ! Set of bodies of the simulation.
    CHARACTER(len=20), ALLOCATABLE, DIMENSION(:) :: BodyFile_list ! List of body files for data input.

    !STEP variables
    INTEGER             :: NSteps     ! Number of time steps of the simulation.
    INTEGER             :: StepCnt    ! Step counter.
    DOUBLE PRECISION    :: TStep      ! Time interval.

    !Aero Matrix Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Full_A_Mat    ! Global Influence Coeficient Matrix.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: Full_RHS       ! Global Right Hand Side Vector.

    !File units/names
    INTEGER           :: FileUnitDB,FileUnitTP,FileUnitRS,FileUnitAERO ! Debug, TecPlot(output), Restart, Results.
    CHARACTER(len=20) :: MVars_FileName ! Main variables input file name

    !Aux/Misc variables
    INTEGER :: ios, ihr, imin, isec, i100th, fhr, fmin, fsec, f100th
    DOUBLE PRECISION :: ti, tf   !Elapsed Time variables

    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U, U2 ! Coordinates (u,theta)
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: V, V2 ! Velocities  (up, thetap)
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE   :: T

    !Aero Force Variables
    DOUBLE PRECISION, DIMENSION(3)    :: Ref_Pt  ! Reference axis point
    DOUBLE PRECISION, DIMENSION(3)    :: Ref_Vrs ! Reference axis versor
	
	!Center of rotation Variables
	DOUBLE PRECISION, DIMENSION(3)    :: axe_Pnt ! Rotation axis point
	DOUBLE PRECISION, DIMENSION(3)    :: axe_Vrs ! Rotation axis versor
	DOUBLE PRECISION, PARAMETER      :: PI=3.1415926535897933 ! Pi.
	INTEGER                          :: auxStpCnt
	DOUBLE PRECISION				 :: p ! T = 2PI/p

    !Main variables input file name
    MVars_FileName='MainV.txt'

    !Main variables read
    call MainVars_Read (MVars_FileName, NSteps, TStep, VinfFile, OPfile, RSfile, DBfile, VinfFlag, OPflag, DBflag, RSflag, RS_time, V_inf, NumberOfBodies, WakeSize, NumberOfGrids, NumberOfBCGrids, NumberOfWakes, BodyFile_list)

    !Output file open
    ExtPnt_OPfile = 'ExtPnt_OPF.dat'
    call TPOP_FileInit (3,OPfile) ! TecPlot output file initialization
    call ExtPts_TPOP_FileInit (8, ExtPnt_OPfile)


    ! Restart
    !Not set yet

    !Bodies construct
    call BodySet_AllocGen (NumberOfBodies, V_inf, BodyFile_list, Body_Set)

    ! Debug
    if (DBflag==1) then
        !Debug file open
        FileUnitDB = 5
        Open (UNIT=FileUnitDB, FILE=DBfile, Status='New', Action='Write', Iostat=ios)
        ! Opening status check
        if (ios /= 0) then
            print * , 'File opening error - Debug File: ',DBfile
            pause
        else
            print * ,'Debug File initialized... '
        end if
    end if

    ! Results File Open (AERO)
    FileUnitAERO=40
    print * ,'pre op AERO'
    pause
    Open (UNIT=FileUnitAERO, FILE='AERO.txt', Status='New', Action='Write', Iostat=ios)
    print * ,'pre op AERO'
    ! Opening status check
    if (ios /= 0) then
        print * , 'File opening error - Result File: AERO.txt'
        pause
    else
        print * ,'Results File "AERO.txt" initialized... '
	end if

    !axe_Pnt  = [-0.225d0, 1.575d0, 0.d0]
	!axe_Pnt  = [-0.225d0, 0.9d0, 0.d0]
 !   axe_Vrs = [0.d0,1.d0,0.d0]
    !
	Ref_Pt  = [0.d0,0.d0,0.d0]
    Ref_Vrs = [0.d0,0.d0,0.d0]
    !
    !print *, 'Moment Reference Axis: Point, Versor ', Ref_Pt , Ref_Vrs
    !print *, 'Vinf', V_inf
    !pause
    
    !call cpu_time(ti)   ! Simulation time initial counter
    CALL GETTIM (ihr, imin, isec, i100th)
    print*,'Initial Time = ',ihr, imin, isec, i100th
    
    DO StepCnt=1,NSteps
        print*,'STEP CNT = ', StepCnt
        !print*, 'grid coords: ', Body_Set(1)%Grids(1)%GridName
		
		!if (StepCnt > 60 .and. StepCnt < 77) then
		!	auxStpCnt = StepCnt - 60
			!if (StepCnt .eq. 66) then
			!	axe_Vrs = [0.d0, 0.d0, 0.d0]
			!end if
		
		      !call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(4),   axe_Pnt, axe_Vrs, 0.d0, 16.d0, 5.d0)
		      !call Coords_Step_Rot (auxStpCnt, TStep, Body_Set(1)%Grids(8),   axe_Pnt, axe_Vrs, 0.d0, 16.d0, 5.d0)
		!end if
		
  !      call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(1),   axe_Pnt, axe_Vrs, 0.d0, 20.d0, 5.d0) !Phase, p, A
		!call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(2),   axe_Pnt, axe_Vrs, 0.d0, 17.d0, 4.d0)
		!call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(3),   axe_Pnt, axe_Vrs, 0.d0, 16.d0, 2.d0)
		!call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(4),   axe_Pnt, axe_Vrs, 0.d0, 19.d0, 3.d0)
		!call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(5),   axe_Pnt, axe_Vrs, 0.d0, 16.d0, 1.d0)
		!call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(6),   axe_Pnt, axe_Vrs, 0.d0, 25.d0, 5.d0)
		!call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(7),   axe_Pnt, axe_Vrs, 0.d0, 18.d0, 3.d0)
		!call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(8),   axe_Pnt, axe_Vrs, 0.d0, 17.d0, 2.d0)
		!
		
		
        !call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%Grids(2),   [0.d0, 20.8d0, 98.2d0], [0.d0 , 523.6d0 , 0.d0])
        !call Coords_Step_Rot (StepCnt, TStep, Body_Set(1)%BCGrids(1), [0.d0, 20.8d0, 98.2d0], [0.d0 , 523.6d0 , 0.d0])

        if (StepCnt == 1) then
            call AeroMat_GenLoc (Body_Set, NumberOfBodies, Full_A_Mat, Full_RHS)
        end if

        call Aero_Step (Body_Set, NumberOfBodies, TStep, StepCnt, NSteps, V_inf, Full_A_Mat, Full_RHS, Ref_Pt, Ref_Vrs)

    END DO
    
    CALL GETTIM (fhr, fmin, fsec, f100th)
    print*,'Initial Time = ',ihr,':', imin,':', isec,'.',i100th
    print*,'Final Time = ',fhr,':', fmin,':', fsec,'.',f100th
    !call cpu_time(tf)   ! Simulation time final counter
    
    tf = dble(fhr)*60.d0*60.d0 + dble(fmin)*60.d0 + dble(fsec) + dble(f100th)*1.d0 /100.d0
    ti = dble(ihr)*60.d0*60.d0 + dble(imin)*60.d0 + dble(isec) + dble(i100th)*1.d0 /100.d0
    print*,'********** ELAPSED TIME = ', tf-ti, '  (sec) **********'

    pause
    end program Main
