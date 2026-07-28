    MODULE STEP
    !Includes all time dependant process as subroutines.

    USE classVortex_Segment
    USE classNodal_Point
    USE classPanel
    USE classGrid
    USE classWake
    USE AeroMatrix
    USE SolverAXB
    USE TecPlot_OutPuts
    USE Aero_Loads
    USE classControl_Point
    USE ShFunct_DShFunct
    USE classBody


    IMPLICIT NONE
    PUBLIC

    CONTAINS ! ==================================================================

    SUBROUTINE Aero_Step (BodySet, nBodies, TStep, StepCnt, NSteps, V_inf, Full_A_Mat, Full_RHS, Ref_Pt, Ref_Vrs)
    ! Defines and controls the aerodynamic step calculation.
    ! Starts with a defined aerodynamic grid and continues with the following flow diagram:
    ! -Generates AeroMatrix and RHS
    ! -Calculates grid panels ring vorticity: Gi. Using the corresponding Solver.
    ! -Produces a step of wake convection. (Induced speed at wake nodes, NewStep, Wake actualization)
    ! -Plots the output wake in default format and TECPLOT format.
    ! The subroutine must be included in a loop for each step required.

    !Input Variables
    INTEGER, INTENT(IN)               :: nBodies ! Number of bodies of the simulation.
    INTEGER, INTENT(IN)               :: StepCnt ! Step counter.
    INTEGER, INTENT(IN)               :: NSteps ! Number Of Steps.
    DOUBLE PRECISION, INTENT(IN)                  :: TStep   ! Time interval.
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: V_inf   ! Free stream velocity (constant or for the current step)
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: Ref_Pt  ! Reference axis point
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: Ref_Vrs ! Reference axis versor

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides structures' data.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:),INTENT(INOUT) :: Full_A_Mat  ! Global Influence Coefficient Matrix
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT (INOUT) :: Full_RHS    ! Global Right Hand Side Vector

    !Local Variables
    DOUBLE PRECISION :: CL, CL1, CL2
    DOUBLE PRECISION :: CD, CD1, CD2
    DOUBLE PRECISION :: CM, CM1, CM2

    INTEGER :: i,j,k,OwAux

    !A Matrix, RHS

    !Ring Vorticities Vector
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: G_vect      ! Global Ring Vorticities Vector

    !Wake adjust ----------------------------------- The step starts adjusting convective nodes before aero matrix calculations.
    print *,'  Wake adjust '
    call ConvectAdjst_Step (BodySet, nBodies, StepCnt)

    !---------------------------------------------------

    !Aero Matrix calculations --------------------------
    !if (StepCnt==1) then
    print *,'A Matrix  '
    call A_Mat_Assembler (BodySet, nBodies, Full_A_Mat)
    !call InfCoefMat_Printing (Full_A_Mat, 40)
    !end if

    print *,' RHS '
    call RHS_Assembler (BodySet, nBodies, V_inf, Full_RHS)
    !call RHSv_Printing (Full_RHS, 40)
    !---------------------------------------------------

    !System Solving ------------------------------------
    call SOLVER (Full_A_Mat, G_vect, Full_RHS)
    !write (FileUnitAERO,*) 'STP n: ',StepCnt
    !write (40,*), 'Step: ', StepCnt
    !call Gv_Printing (G_vect, 40)
    !---------------------------------------------------

    !Vorticity processing-------------------------------
    print *,'  Vorticity processing '
    call Vorticity_Process (BodySet, nBodies, G_vect)
    !---------------------------------------------------

    !Induced Speed calculations ------------------------
    print *,'  InducedSpeed Step '
    call InducedSpeed_Step (BodySet, nBodies)

    !if (StepCnt == NSteps) then
    !print *,'  ExtraPoint processing '
    !call IndSpeed_atExtPoint (BodySet(1)%ExtraPoints, BodySet, nBodies)
    !call ExtPt_Vel (BodySet(1)%ExtraPoints, V_inf)
    !call ExtPnt_TPOP (BodySet(1)%ExtraPoints, StepCnt, StepCnt, 8)
    !end if


    !---------------------------------------------------

    !Aerodynamic Loads----------------------------------
    print *,'  AeroLoads '
    call DeltaCp_calc (BodySet, nBodies, TStep, V_inf, StepCnt)
    call ForceCoeficient_calc (BodySet(1)%Grids(1), V_inf, CL, CD, CM, Ref_Pt, Ref_Vrs)
    !call ForceCoeficient_calc (BodySet(1)%Grids(2), V_inf, CL1, CD1, CM1, Ref_Pt, Ref_Vrs)
    !call ForceCoeficient_calc (BodySet(1)%Grids(3), V_inf, CL1, CD1)
    !call ForceCoeficient_calc (BodySet(2)%Grids(1), V_inf, CL2, CD2)
    !write (20,*) 'STP n: ',StepCnt
    !write (40,'(I3,A2,F10.7,A2,F10.7,A2,F10.7,A2,F10.7,A2,F10.7,A2,F10.7)'), StepCnt,': ', CL,'  ', CD,'  ', CL1,'  ', CD1,'  ', CM,'  ', CM1
    write (40,'(I3,A2,F10.7,A2,F10.7,A2,F10.7)'), StepCnt,': ', CL,'  ', CD,'  ', CM
    !write (5,'(A8,I4)'), 'Step: ', StepCnt
	!if (StepCnt == NSteps) then
		!do i=1,8
		!	call PanelSet_Printing (BodySet(1)%Grids(i)%PanelSet, 40)
		!end do
		!write (40,*)
	!end if
	
	!call GridPrinting (BodySet(1)%Grids(1), 40)
    !write (5,'(I3,A2,F10.7,A2,F10.7,A2,F10.7,A2,F10.7)'), StepCnt,': ', CL,'  ', CD,'  ',CL2,'  ', CD2
    !print *, 'DeltaCP check _: ', StepCnt, BodySet(1)%Grids(1)%PanelSet(1)%CtrlPoint%Delta_Cp
    !---------------------------------------------------

    !Wake Convection -----------------------------------
    call CONVECTION (BodySet, nBodies, TStep, StepCnt, V_inf)
    !---------------------------------------------------

    !TPOP ----------------------------------------------
    !if (StepCnt == NSteps) then
        call TPOP_Step (BodySet,nBodies,StepCnt,3)
    !end if
    !---------------------------------------------------
    print *,'  End Step '
    END SUBROUTINE Aero_Step

    SUBROUTINE Coords_Step_Rot (StepCnt, TStep, GRD, axe_Pnt, axe_Vrs, f, p, A_deg)!Phase, p, A
    !Gives the grid (or BC_grid) nodes a prescribed rotation. In addition, gives the Control Points its own velocity associated.
    !Actualices the nodes coordinates.

    interface
    FUNCTION cross_product(U,V)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
    DOUBLE PRECISION, DIMENSION(3)              :: cross_product
    DOUBLE PRECISION                            :: cross1, cross2, cross3

    END FUNCTION cross_product
    end interface

    !Input Variables
    INTEGER, INTENT(IN) :: StepCnt
    DOUBLE PRECISION, INTENT(IN)    :: TStep
	DOUBLE PRECISION, INTENT(IN)    :: f  ! angle face
	DOUBLE PRECISION, INTENT(IN)    :: A_deg  ! angle Amplitud (deg)
	DOUBLE PRECISION, INTENT(IN)    :: p  ! Period: T= 2PI/p
	
    !DOUBLE PRECISION, DIMENSION (3), INTENT(IN)    :: rot_ctr !Center of rotation
    !DOUBLE PRECISION, DIMENSION (3), INTENT(INOUT) :: OM ! Vector Angular Velocity
	DOUBLE PRECISION, DIMENSION (3), INTENT(IN)    :: axe_Pnt !Center of rotation
    DOUBLE PRECISION, DIMENSION (3), INTENT(IN)    :: axe_Vrs ! Vector Angular Velocity

    !Output Variables
    CLASS(GRID), INTENT(INOUT)   :: GRD !Grid to impose the rotation - POLYMOROHIC VARIABLE

    !Local Variables
    !DOUBLE PRECISION, DIMENSION (3) :: OM ! Vector Angular Velocity
    DOUBLE PRECISION, DIMENSION (3) :: Vv ! Vector Linear Velocity
    DOUBLE PRECISION, DIMENSION (3) :: X_temp, R, A! Temp Coordinates vector
    DOUBLE PRECISION, DIMENSION (3) :: OwnV ! Own velocity of a control point
    DOUBLE PRECISION                :: r_abs0, r_abs, th_0
    INTEGER  :: nNodes ! Number of nodes of the grid
    INTEGER  :: nPnls ! Number of panels of the grid
    INTEGER  :: i,j      ! Loop index
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Xa
    DOUBLE PRECISION, DIMENSION (8)   :: PShapeFunct
    DOUBLE PRECISION, DIMENSION (8,2) :: Derived_PSF
    DOUBLE PRECISION, DIMENSION (3)   :: NrmlVersor, XIv, CP_xyz
	DOUBLE PRECISION, DIMENSION (3)  :: OM
	DOUBLE PRECISION  :: beta
	DOUBLE PRECISION  :: Amp  ! angle Amplitud (rad)
	DOUBLE PRECISION, PARAMETER      :: PI=3.1415926535897933 ! Pi.

    !beta = 0.087266462599716d0*sin(2*PI*StepCnt/16.d0 + f) ! \beta, rotation angle
	!OM =  axe_Vrs * 0.087266462599716d0*(2*PI/16.d0)*cos(2*PI*StepCnt/16.d0 + f) ! \omega, angular velocity
	
	Amp = A_deg*PI/180.d0; ! deg2rad()
	
	beta = Amp*sin(2*PI*StepCnt/p + f) ! \beta, rotation angle
	OM =  axe_Vrs * Amp*(2*PI/p)*cos(2*PI*StepCnt/p + f) ! \omega, angular velocity

    nNodes = GRD%nNodalPoints
    nPnls =  GRD%nPanels

    !Nodes Position Update
    DO i=1,nNodes
		if (GRD%NodalPointSet(i)%xyz(1) < axe_Pnt(1)) then
		
			X_temp = GRD%NodalPointSet(i)%xyz - axe_Pnt ! pt
			r_abs = norm2(dot_product(X_temp,axe_Vrs)*axe_Vrs - X_temp)
		
			A(1) = (1.d0 - cos(beta))*r_abs
			A(2) = 0.d0
			A(3) = sin(beta)*r_abs
		
			GRD%NodalPointSet(i)%xyz = GRD%NodalPointSet(i)%xyz_loc + A
		end if

        !GRD%NodalPointSet(i)%xyz(1) = (- OM(3) * X_temp(2) * TStep + GRD%NodalPointSet(i)%xyz(1)) !* (norm2(GRD%NodalPointSet(i)%xyz_loc) / norm2(GRD%NodalPointSet(i)%xyz))
        !GRD%NodalPointSet(i)%xyz(2) = (  OM(3) * X_temp(1) * TStep + GRD%NodalPointSet(i)%xyz(2)) !* (norm2(GRD%NodalPointSet(i)%xyz_loc) / norm2(GRD%NodalPointSet(i)%xyz))
        !GRD%NodalPointSet(i)%xyz(1) = GRD%NodalPointSet(i)%xyz(1) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))
        !GRD%NodalPointSet(i)%xyz(2) = GRD%NodalPointSet(i)%xyz(2) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))

        !A = cross_product (OM, X_temp) * TStep

        !A(1) = - OM(3) * X_temp(2) * TStep
        !A(2) = OM(3) * X_temp(1) * TStep
        !A(3) = 0.0

        !GRD%NodalPointSet(i)%xyz(1) = ( A(1) + X_temp(1)) * (norm2(X_temp) / norm2(A + X_temp)) + rot_ctr(1)
        !GRD%NodalPointSet(i)%xyz(2) = ( A(2) + X_temp(2))  * (norm2(X_temp) / norm2(A + X_temp)) + rot_ctr(2)

        !GRD%NodalPointSet(i)%xyz = (A + X_temp) * (norm2(X_temp) / norm2(A + X_temp)) + rot_ctr

        !GRD%NodalPointSet(i)%xyz(1) = GRD%NodalPointSet(i)%xyz(1) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))
        !GRD%NodalPointSet(i)%xyz(2) = GRD%NodalPointSet(i)%xyz(2) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))


        !end if
        !print *, 'Nc1',GRD%NodalPointSet(i)%xyz(1)
        !print *, 'Nc2',GRD%NodalPointSet(i)%xyz(2)
        !print *, 'Nc3',GRD%NodalPointSet(i)%xyz(3)
        !pause
    END DO

    !CPs Speed
    DO j=1,nPnls
		if (GRD%PanelSet(j)%CtrlPoint%xyz(1) < axe_Pnt(1)) then
		
			X_temp = GRD%PanelSet(j)%CtrlPoint%xyz - axe_Pnt ! pt
			R = dot_product(X_temp,axe_Vrs)*axe_Vrs - X_temp !r vector

			!OwnV(1) = - OM(3) * (GRD%PanelSet(j)%CtrlPoint%xyz(2) - rot_ctr (2))
			!OwnV(2) =   OM(3) * (GRD%PanelSet(j)%CtrlPoint%xyz(1) - rot_ctr (1))
			!OwnV(3) =   0.0

			OwnV =  cross_product (OM, R)

			!print *, 'pnl: ',j
			!print *, 'OV1',OwnV(1)
			!print *, 'OV2',OwnV(2)
			!print *, 'OV3',OwnV(3)
			!pause
			call CP_Velocity (GRD%PanelSet(j)%CtrlPoint, GRD%PanelSet(j)%CtrlPoint%WindVelocity, OwnV)
		end if
    END DO

    !CPs Processing
    DO j=1,nPnls
        call InGrid_Xa (GRD%PanelSet(j)%PanelNodes, GRD%NodalPointSet, GRD%PanelSet(j)%Panel_Type, Xa)
        XIv = [0.0,0.0,0.0]
        !print *, 'Xa ',Xa
        !pause
        call Shape_Functions (GRD%PanelSet(j)%Panel_Type, PShapeFunct, XIv)
        call Derived_Shape_Functions (GRD%PanelSet(j)%Panel_Type, Derived_PSF)
        call CP_Coords (Xa, PShapeFunct, CP_xyz)
        call CP_NormalVersor (Xa, Derived_PSF, NrmlVersor)
        !print *, 'CpCorr, Nvrsr',CP_xyz, NrmlVersor
        !    pause
        GRD%PanelSet(j)%CtrlPoint%xyz=CP_xyz
        GRD%PanelSet(j)%CtrlPoint%NormalVersor=NrmlVersor
    END DO

	END SUBROUTINE Coords_Step_Rot
	
	!SUBROUTINE Coords_Step_Rot (StepCnt, TStep, GRD, rot_ctr, OM)  ------------->>>ORIGINAL
 !   !Gives the grid (or BC_grid) nodes a prescribed rotation. In addition, gives the Control Points its own velocity associated.
 !   !Actualices the nodes coordinates.
 !
 !   interface
 !   FUNCTION cross_product(U,V)
 !
 !   IMPLICIT NONE
 !
 !   DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
 !   DOUBLE PRECISION, DIMENSION(3)              :: cross_product
 !   DOUBLE PRECISION                            :: cross1, cross2, cross3
 !
 !   END FUNCTION cross_product
 !   end interface
 !
 !   !Input Variables
 !   INTEGER, INTENT(IN) :: StepCnt
 !   DOUBLE PRECISION, INTENT(IN)    :: TStep
 !   DOUBLE PRECISION, DIMENSION (3), INTENT(IN)    :: rot_ctr !Center of rotation
 !   DOUBLE PRECISION, DIMENSION (3), INTENT(INOUT) :: OM ! Vector Angular Velocity
 !
 !   !Output Variables
 !   CLASS(GRID), INTENT(INOUT)   :: GRD !Grid to impose the rotation - POLYMOROHIC VARIABLE
 !
 !   !Local Variables
 !   !DOUBLE PRECISION, DIMENSION (3) :: OM ! Vector Angular Velocity
 !   DOUBLE PRECISION, DIMENSION (3) :: Vv ! Vector Linear Velocity
 !   DOUBLE PRECISION, DIMENSION (3) :: X_temp, A, B! Temp Coordinates vector
 !   DOUBLE PRECISION, DIMENSION (3) :: OwnV ! Own velocity of a control point
 !   DOUBLE PRECISION                :: r_abs0, r_abs, th_0
 !   INTEGER  :: nNodes ! Number of nodes of the grid
 !   INTEGER  :: nPnls ! Number of panels of the grid
 !   INTEGER  :: i,j      ! Loop index
 !   DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Xa
 !   DOUBLE PRECISION, DIMENSION (8)   :: PShapeFunct
 !   DOUBLE PRECISION, DIMENSION (8,2) :: Derived_PSF
 !   DOUBLE PRECISION, DIMENSION (3)   :: NrmlVersor, XIv, CP_xyz
 !
 !   !OM = [0.0 , 0.0 , -20.0]
 !
 !   nNodes = GRD%nNodalPoints
 !   nPnls =  GRD%nPanels
 !
 !   !Nodes Position Update
 !   DO i=1,nNodes
	!	if GRD%NodalPointSet(i)%xyz(1) < rot_ctr(1) then
	!	
 !       X_temp = GRD%NodalPointSet(i)%xyz - rot_ctr
 !
 !       !GRD%NodalPointSet(i)%xyz(1) = (- OM(3) * X_temp(2) * TStep + GRD%NodalPointSet(i)%xyz(1)) !* (norm2(GRD%NodalPointSet(i)%xyz_loc) / norm2(GRD%NodalPointSet(i)%xyz))
 !       !GRD%NodalPointSet(i)%xyz(2) = (  OM(3) * X_temp(1) * TStep + GRD%NodalPointSet(i)%xyz(2)) !* (norm2(GRD%NodalPointSet(i)%xyz_loc) / norm2(GRD%NodalPointSet(i)%xyz))
 !       !GRD%NodalPointSet(i)%xyz(1) = GRD%NodalPointSet(i)%xyz(1) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))
 !       !GRD%NodalPointSet(i)%xyz(2) = GRD%NodalPointSet(i)%xyz(2) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))
 !
 !       A = cross_product (OM, X_temp) * TStep
 !
 !       !A(1) = - OM(3) * X_temp(2) * TStep
 !       !A(2) = OM(3) * X_temp(1) * TStep
 !       !A(3) = 0.0
 !
 !       !GRD%NodalPointSet(i)%xyz(1) = ( A(1) + X_temp(1)) * (norm2(X_temp) / norm2(A + X_temp)) + rot_ctr(1)
 !       !GRD%NodalPointSet(i)%xyz(2) = ( A(2) + X_temp(2))  * (norm2(X_temp) / norm2(A + X_temp)) + rot_ctr(2)
 !
 !       GRD%NodalPointSet(i)%xyz = (A + X_temp) * (norm2(X_temp) / norm2(A + X_temp)) + rot_ctr
 !
 !       !GRD%NodalPointSet(i)%xyz(1) = GRD%NodalPointSet(i)%xyz(1) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))
 !       !GRD%NodalPointSet(i)%xyz(2) = GRD%NodalPointSet(i)%xyz(2) * (norm2(GRD%NodalPointSet(i)%xyz_loc - rot_ctr) / norm2(GRD%NodalPointSet(i)%xyz - rot_ctr))
 !
 !
 !       !end if
 !       !print *, 'Nc1',GRD%NodalPointSet(i)%xyz(1)
 !       !print *, 'Nc2',GRD%NodalPointSet(i)%xyz(2)
 !       !print *, 'Nc3',GRD%NodalPointSet(i)%xyz(3)
 !       !pause
 !   END DO
 !
 !   !CPs Speed
 !   DO j=1,nPnls
 !
 !       !OwnV(1) = - OM(3) * (GRD%PanelSet(j)%CtrlPoint%xyz(2) - rot_ctr (2))
 !       !OwnV(2) =   OM(3) * (GRD%PanelSet(j)%CtrlPoint%xyz(1) - rot_ctr (1))
 !       !OwnV(3) =   0.0
 !
 !       OwnV =  cross_product (OM, (GRD%PanelSet(j)%CtrlPoint%xyz - rot_ctr))
 !
 !       !print *, 'pnl: ',j
 !       !print *, 'OV1',OwnV(1)
 !       !print *, 'OV2',OwnV(2)
 !       !print *, 'OV3',OwnV(3)
 !       !pause
 !       call CP_Velocity (GRD%PanelSet(j)%CtrlPoint, GRD%PanelSet(j)%CtrlPoint%WindVelocity, OwnV)
 !   END DO
 !
 !   !CPs Processing
 !   DO j=1,nPnls
 !       call InGrid_Xa (GRD%PanelSet(j)%PanelNodes, GRD%NodalPointSet, GRD%PanelSet(j)%Panel_Type, Xa)
 !       XIv = [0.0,0.0,0.0]
 !       !print *, 'Xa ',Xa
 !       !pause
 !       call Shape_Functions (GRD%PanelSet(j)%Panel_Type, PShapeFunct, XIv)
 !       call Derived_Shape_Functions (GRD%PanelSet(j)%Panel_Type, Derived_PSF)
 !       call CP_Coords (Xa, PShapeFunct, CP_xyz)
 !       call CP_NormalVersor (Xa, Derived_PSF, NrmlVersor)
 !       !print *, 'CpCorr, Nvrsr',CP_xyz, NrmlVersor
 !       !    pause
 !       GRD%PanelSet(j)%CtrlPoint%xyz=CP_xyz
 !       GRD%PanelSet(j)%CtrlPoint%NormalVersor=NrmlVersor
 !   END DO
 !
 !   END SUBROUTINE Coords_Step_Rot

    SUBROUTINE Coords_Step (StepCnt, TStep, GRD, U, V, T)
    !Gives the grid nodes a prescribed time variable motion in each step. In addition, gives the Control Points its own velocity associated.
    !Actualices the nodes coordinates.

    !Input Variables
    INTEGER, INTENT(IN) :: StepCnt
    DOUBLE PRECISION, INTENT(IN)    :: TStep
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE, INTENT(IN) :: U ! Coordinates (u,theta)
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE, INTENT(IN) :: V ! Velocities  (up, thetap)
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE, INTENT(IN)   :: T ! Dimensional time steps

    !Output Variables
    type(Grid),INTENT(INOUT)    :: GRD ! Input Grid

    !Local Variables
    !DOUBLE PRECISION     :: OmU ! Plunge frequency
    !DOUBLE PRECISION     :: AU  ! Plunge amplitude
    !DOUBLE PRECISION     :: OmT ! Pitch frequency
    !DOUBLE PRECISION     :: AT  ! Pitch amplitude
    DOUBLE PRECISION     :: l, lx, ly, phi
    DOUBLE PRECISION, DIMENSION (2) :: Ax
    DOUBLE PRECISION, DIMENSION(3) :: OwnV ! Own velocity of a control point
    INTEGER  :: nNodes ! Number of nodes of the grid
    INTEGER  :: nPnls ! Number of panels of the grid
    !DOUBLE PRECISION     :: dU     ! Plunge Coodinate
    !DOUBLE PRECISION     :: dT     ! Pitch Coordinate
    !DOUBLE PRECISION     :: Up     ! Plunge speed
    !DOUBLE PRECISION     :: Tp     ! Pitch speed
    DOUBLE PRECISION     :: c      ! Chord
    INTEGER  :: i,j      ! Loop index
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Xa
    DOUBLE PRECISION, DIMENSION (8)   :: PShapeFunct
    DOUBLE PRECISION, DIMENSION (8,2) :: Derived_PSF
    DOUBLE PRECISION, DIMENSION (3)   :: NrmlVersor, XIv, CP_xyz

    c=6.0

    nNodes = GRD%nNodalPoints
    nPnls =  GRD%nPanels



    !DO i=1,nNodes
    !    GRD%NodalPointSet(i)%xyz(1) = (GRD%NodalPointSet(i)%xyz_loc(1))*cos(U(StepCnt,2))!+c/2
    !    GRD%NodalPointSet(i)%xyz(3) = GRD%NodalPointSet(i)%xyz(3) - (GRD%NodalPointSet(i)%xyz_loc(1))*sin(U(StepCnt,2))
    !END DO

    !DO i=1,nNodes
    !    lx=sqrt((GRD%NodalPointSet(i)%xyz_loc(3))**2 + (GRD%NodalPointSet(i)%xyz_loc(2))**2)!*sin(U(StepCnt,2)*0.5)
    !
    !    GRD%NodalPointSet(i)%xyz(2) = - lx * cos(U(StepCnt,1))
    !    !GRD%NodalPointSet(i)%xyz(3) = lx * sin(U(StepCnt,2))
    !
    !    ly=sqrt((GRD%NodalPointSet(i)%xyz_loc(1))**2 + (GRD%NodalPointSet(i)%xyz_loc(3))**2)!*sin(U(StepCnt,2)*0.5)
    !
    !    GRD%NodalPointSet(i)%xyz(1) = ly * sin(U(StepCnt,2))
    !    GRD%NodalPointSet(i)%xyz(3) = lx * sin(U(StepCnt,1)) + ly * cos(U(StepCnt,2))
    !
    !    !print *, 'Nc1',GRD%NodalPointSet(i)%xyz(1)
    !    !print *, 'Nc2',GRD%NodalPointSet(i)%xyz(2)
    !    !print *, 'Nc3',GRD%NodalPointSet(i)%xyz(3)
    !    !pause
    !END DO

    Ax=[2.0 , 53.0]
    DO i=1,nNodes
        l=sqrt((GRD%NodalPointSet(i)%xyz_loc(1)-Ax(1))**2 + (GRD%NodalPointSet(i)%xyz_loc(2)-Ax(2))**2)!*sin(U(StepCnt,2)*0.5)
        if (GRD%NodalPointSet(i)%xyz_loc(1)-Ax(1) == 0.0) then
            phi = 0.0
        else
            phi = atan((GRD%NodalPointSet(i)%xyz_loc(2)-Ax(2))/(GRD%NodalPointSet(i)%xyz_loc(1)-Ax(1)))+(U(StepCnt,2))
        end if
        !print *, 'phi',phi
        !GRD%NodalPointSet(i)%xyz(1) = GRD%NodalPointSet(i)%xyz(1) + l*sin(U(StepCnt,2))!+c/2
        !GRD%NodalPointSet(i)%xyz(2) = GRD%NodalPointSet(i)%xyz(2) + l*cos(U(StepCnt,2))
        GRD%NodalPointSet(i)%xyz(1) = Ax(1) + l * (cos(phi)) * sign(1.0 , GRD%NodalPointSet(i)%xyz_loc(1) - Ax(1))
        GRD%NodalPointSet(i)%xyz(2) = Ax(2) + l * (sin(phi)) * sign(1.0 , GRD%NodalPointSet(i)%xyz_loc(1) - Ax(1))
        !print *, 'Nc1',GRD%NodalPointSet(i)%xyz(1)
        !print *, 'Nc2',GRD%NodalPointSet(i)%xyz(2)
        !print *, 'Nc3',GRD%NodalPointSet(i)%xyz(3)
        !pause
    END DO


    !CPs Speed
    DO j=1,nPnls
        !if (GRD%PanelSet(j)%CtrlPoint%xyz_loc(1)-Ax(1) == 0.0) then
        !   phi = 0.0
        !else
        !    phi = atan((GRD%PanelSet(j)%CtrlPoint%xyz_loc(2)-Ax(2))/(GRD%PanelSet(j)%CtrlPoint%xyz_loc(1)-Ax(1)))+(U(StepCnt,2))
        !end if
        !
        !OwnV(1) = sqrt((GRD%PanelSet(j)%CtrlPoint%xyz_loc(1)-Ax(1))**2 + (GRD%PanelSet(j)%CtrlPoint%xyz_loc(2)-Ax(2))**2)*V(StepCnt,2)*sin(phi) * sign(1.0 , Ax(1) - GRD%PanelSet(j)%CtrlPoint%xyz_loc(1))
        !OwnV(2) = sqrt((GRD%PanelSet(j)%CtrlPoint%xyz_loc(1)-Ax(1))**2 + (GRD%PanelSet(j)%CtrlPoint%xyz_loc(2)-Ax(2))**2)*V(StepCnt,2)*cos(phi) * sign(1.0 , Ax(2) - GRD%PanelSet(j)%CtrlPoint%xyz_loc(2))
        !OwnV(3) = 0.0

        !OwnV(1) = - (GRD%PanelSet(j)%CtrlPoint%xyz_loc(1))*V(StepCnt,2)*sin(U(StepCnt,2))
        !OwnV(2) = 0.0
        !OwnV(3) = V(StepCnt,1) - (GRD%PanelSet(j)%CtrlPoint%xyz_loc(1))*V(StepCnt,2)*cos(U(StepCnt,2))

        lx=sqrt((GRD%PanelSet(j)%CtrlPoint%xyz_loc(3))**2 + (GRD%PanelSet(j)%CtrlPoint%xyz_loc(2))**2)
        ly=sqrt((GRD%PanelSet(j)%CtrlPoint%xyz_loc(3))**2 + (GRD%PanelSet(j)%CtrlPoint%xyz_loc(1))**2)

        OwnV(1) = - V(StepCnt,2) * ly * sin(U(StepCnt,2))
        OwnV(2) = V(StepCnt,1) * lx * sin(U(StepCnt,1))
        OwnV(3) = V(StepCnt,1) * lx * cos(U(StepCnt,1)) + V(StepCnt,2) * ly * cos(U(StepCnt,2))

        !print *, 'pnl: ',j
        !print *, 'OV1',OwnV(1)
        !print *, 'OV2',OwnV(2)
        !print *, 'OV3',OwnV(3)

        call CP_Velocity (GRD%PanelSet(j)%CtrlPoint, GRD%PanelSet(j)%CtrlPoint%WindVelocity, OwnV)
    END DO

    !CPs Processing
    DO j=1,nPnls
        call InGrid_Xa (GRD%PanelSet(j)%PanelNodes, GRD%NodalPointSet, GRD%PanelSet(j)%Panel_Type, Xa)
        XIv = [0.0,0.0,0.0]
        !print *, 'Xa ',Xa
        !pause
        call Shape_Functions (GRD%PanelSet(j)%Panel_Type, PShapeFunct, XIv)
        call Derived_Shape_Functions (GRD%PanelSet(j)%Panel_Type, Derived_PSF)
        call CP_Coords (Xa, PShapeFunct, CP_xyz)
        call CP_NormalVersor (Xa, Derived_PSF, NrmlVersor)
        !print *, 'CpCorr, Nvrsr',CP_xyz, NrmlVersor
        !    pause
        GRD%PanelSet(j)%CtrlPoint%xyz=CP_xyz
        GRD%PanelSet(j)%CtrlPoint%NormalVersor=NrmlVersor
    END DO

    END SUBROUTINE Coords_Step


    SUBROUTINE Motion_Input (U,V,T,TStep)
    !Subroutine for loading the grid motion as arrays of temporal variables for plunge and pitch coordinates.

    !Input Variables
    DOUBLE PRECISION, INTENT(IN) :: TStep

    !Output Variables
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE, INTENT(INOUT) :: U ! Coordinates (u,theta)
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE, INTENT(INOUT) :: V ! Velocities  (up, thetap)
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE, INTENT(INOUT)   :: T ! Dimensional time steps
    !Local Variables
    CHARACTER (len=40)  :: FileName
    INTEGER             :: Steps
    INTEGER             :: i, ios
    DOUBLE PRECISION, DIMENSION(5)  :: Line

    !Number of Time steps
    Steps = 200

    !Array allocation
    allocate(U(Steps,2))
    allocate(V(Steps,2))
    allocate(T(Steps))

    !File open
    FileName = 'DRS_mot.dat'
    !FileName = 'Mot_input_fly.dat'
    !FileName = 'desp_vel_full.txt'
    open (1, FILE=FileName, IOSTAT=ios, ACTION='read')

    !DRS, Fly
    DO i=1,Steps
        read (1,*,IOSTAT=ios) Line
        T(i)   = Line(1)
        U(i,1) = Line(2)  !u
        U(i,2) = Line(3)  !theta
        V(i,1) = Line(4)  !u_p
        V(i,2) = Line(5)  !theta_p
    END DO

    !Marcos
    !DO i=1,Steps
    !    read (1,*,IOSTAT=ios) Line
    !    T(i)   = Line(1)
    !    U(i,1) = -Line(2)/10.0                      !u
    !    U(i,2) = Line(3)                       !theta  0.0637351
    !    V(i,1) = -Line(4)/156.9            !u_p
    !    V(i,2) = Line(5)*10.0/156.9             !theta_p
    !END DO

    print *, 'T : ',T
    print *, 'U : ',U
    print *, 'V : ',V

    close (UNIT=1)



    END SUBROUTINE Motion_Input


    SUBROUTINE CONVECTION (BodySet, nBodies, TStep, StepCnt, V_inf)
    ! Defines and organices the convection process for a Grid with ConvectLine into a Wake and goes through every wake in the BodySet.
    ! WkeSize and TStep are the number of convection steps that will be stored and their time step, respectivly.
    ! Includes the specific routines in this module, asuming that the wake is already allocatad and initiallizated.
    ! Assumes that the Wake and its Owners (ConvectLine and Grid) are from the same body.
    ! Requires previous Induced Speed at Wake Nodes calculation.

    !Input Variables
    INTEGER, INTENT(IN)               :: nBodies
    INTEGER, INTENT(IN)               :: StepCnt
    DOUBLE PRECISION, INTENT(IN)                  :: TStep
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: V_inf

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.

    !Local Variables
    type(Nodal_Point), ALLOCATABLE, DIMENSION(:)    :: NewNodes ! Group of nodes generated in each step of convection, forming a "parallel" line from convection line.
    type(Vortex_Segment), ALLOCATABLE, DIMENSION(:) :: NewSegments !Group of segments genereated according to the 'NewNodes' array.
    INTEGER     :: OwAux  ! Owner Grid extraction variable
    INTEGER     :: i,j,k  ! Loop Index


    DO i=1,nBodies
        DO j=1,BodySet(i)%nWakes

            WKEj: associate (WKE => BodySet(i)%Wakes(j))

                !Owner Search
                DO k=1, BodySet(i)%nGrids
                    if (BodySet(i)%Grids(k)%GridLabel == WKE%Owner) then
                        OwAux = k
                        print *, 'OwAux  ', k
                        exit
                    end if
                END DO


                !call ConvectNodes_Adjust (WKE, BodySet(i)%Grids(OwAux), WKE%ConvecLine, StepCnt)

                call Wake_NewStep (WKE%ConvecLine, WKE%WakeLabel, BodySet(i)%Grids(OwAux), TStep, NewNodes, NewSegments, V_inf)   ! Generates generic new nodes and segments as arrays.

                call Wake_Actualization (WKE, NewNodes, NewSegments, StepCnt, WKE%WakeSize, TStep, V_inf) ! Inserts the new nodes and segments in the existing wake and actualices the position and indexes of the previous ones.

            end associate WKEj

        END DO
    END DO

    END SUBROUTINE CONVECTION

    SUBROUTINE ConvectAdjst_Step (BodySet, nBodies, StepCnt)
    ! Actualices the position of the last loaded nodes of wake according to the grid movement.

    !Input Variables
    INTEGER, INTENT(IN)               :: nBodies
    INTEGER, INTENT(IN)               :: StepCnt

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.

    !Local Variables
    INTEGER     :: OwAux  ! Owner Grid extraction variable
    INTEGER     :: i,j,k  ! Loop Index


    DO i=1,nBodies
        DO j=1,BodySet(i)%nWakes

            WKEj: associate (WKE => BodySet(i)%Wakes(j))

                !Owner Search
                DO k=1, BodySet(i)%nGrids
                    if (BodySet(i)%Grids(k)%GridLabel == WKE%Owner) then
                        OwAux = k
                        print *, 'OwAux  ', k
                        exit
                    end if
                END DO


                call ConvectNodes_Adjust (WKE, BodySet(i)%Grids(OwAux), WKE%ConvecLine, StepCnt)

            end associate WKEj

        END DO
    END DO






    END SUBROUTINE ConvectAdjst_Step

    END MODULE STEP
