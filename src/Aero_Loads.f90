    MODULE Aero_Loads

    !Gathers subroutines for Aerodynamics Loads Computations

    USE classVortex_Segment
    USE classNodal_Point
    USE classPanel
    USE classPoint
    USE classControl_Point
    USE classGrid
    USE classWake
    USE classConvection_Line
    USE classBC_Grid
    USE classBody

    IMPLICIT NONE
    PUBLIC

    CONTAINS

    SUBROUTINE InducedSpeed_Step (BodySet, nBodies)
    ! Calculates the induced speed at every Point (Nodes, CtrlPoints) for further Aero_Loads calculations.
    ! Executes once at every step and loads the results into the Points' IndSpeed variables.

    !Input Variables
    INTEGER, INTENT(IN) :: nBodies

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: BodySet

    !Local Variables
    INTEGER :: i,j, k, l !Loop index

    DO i=1,nBodies !Bodies loop
        DO j=1,BodySet(i)%nGrids !Grids loop
            !print *, 'Indspeed grids CP---'
            call IndSpeed_atCtrlPoints (BodySet(i)%Grids(j), BodySet, nBodies) !Calculates indspeed at Grids' CtrlPoints
            !pause
        END DO
        DO k=1,BodySet(i)%nBCGrids !BCGrids loop
            !print *, 'Indspeed BCgrids CP---'
            call IndSpeed_atCtrlPoints (BodySet(i)%BCGrids(k), BodySet, nBodies) !Calculates indspeed at BCGrids' CtrlPoints
            !pause
        END DO
        DO l=1,BodySet(i)%nWakes !Wakes loop
            !print *, 'Indspeed Wakes---'
            call Wake_NodesIndSpeed (BodySet(i)%Wakes(l), BodySet, nBodies) !Calculates indspeed at Wakes' NodalPoints
            !pause
        END DO
    END DO


    END SUBROUTINE InducedSpeed_Step

    SUBROUTINE Wake_NodesIndSpeed (WKE, BodySet, nBodies)
    !Calculates the induced speed at every Node of the Wake (WKE).
    !Considering all Wakes and Grids (BC_Grids as well) influence.
    !Free stream speed (Vinf) is NOT considered.

    !Input Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    INTEGER, INTENT(IN)                                 :: nBodies ! Number of bodies of the simulation.

    !Output Variables
    type(WAKE), INTENT(INOUT)   :: WKE   ! Provides the nodes for induced speed calculation.


    !Loca Variables
    INTEGER :: nNodes
    INTEGER :: i

    nNodes = WKE%nNodalPoints

    !DO i=1,nNodes
    if (WKE%LLStride < 0 .or. WKE%LLStride==WKE%WakeSize-1) then
    
        !$OMP PARALLEL SHARED (WKE, BodySet, nBodies)
        !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(i)
    
        DO i=1,WKE%NodesStride*abs(WKE%LLStride+1)+WKE%NodesStride
            !if (WKE%NodalPointSet(i)%Label == 0) then
            !exit
            !else
            call IndSpeed_atNode (BodySet, nBodies, WKE%NodalPointSet(i))
            !end if
        END DO
        
        !$OMP END DO 
        !$OMP END PARALLEL
    else
    
       !$OMP PARALLEL SHARED (WKE, BodySet, nBodies)
       !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(i)
        DO i=1,nNodes
            call IndSpeed_atNode (BodySet, nBodies, WKE%NodalPointSet(i))
        END DO
        
       !$OMP END DO 
       !$OMP END PARALLEL
        
    end if


    END SUBROUTINE Wake_NodesIndSpeed
    
    !SUBROUTINE Wake_NodesIndSpeed_old (WKE, BodySet, nBodies)
    !!Calculates the induced speed at every Node of the Wake (WKE).
    !!Considering all Wakes and Grids (BC_Grids as well) influence.
    !!Free stream speed (Vinf) is NOT considered.
    !
    !!Input Variables
    !type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    !INTEGER, INTENT(IN)                                 :: nBodies ! Number of bodies of the simulation.
    !
    !!Output Variables
    !type(WAKE), INTENT(INOUT)   :: WKE   ! Provides the nodes for induced speed calculation.
    !
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
    !        call IndSpeed_atNode (BodySet, nBodies, WKE%NodalPointSet(i))
    !        !end if
    !    END DO
    !else
    !    DO i=1,nNodes
    !        call IndSpeed_atNode (BodySet, nBodies, WKE%NodalPointSet(i))
    !    END DO
    !end if
    !
    !
    !END SUBROUTINE Wake_NodesIndSpeed_old

    SUBROUTINE IndSpeed_atCtrlPoints (GRD, BodySet, nBodies)

    !Calculates the induced speed at every Control Point of the Grid (GRD).
    !Considering all Wakes and Grids (BC_Grids as well) influence.
    !Free stream speed (Vinf) is NOT considered.
    ! USES POLYMORPHIC VARIABLES SO AS TO PROCESS GRIDS AND BC_GRIDS.

    !Input Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    INTEGER, INTENT(IN)                                 :: nBodies ! Number of bodies of the simulation.

    !Output Variables
    CLASS(GRID), INTENT(INOUT)   :: GRD   !Provides the control points for induced speed calculation. POLYMORPHIC FOR GIRDS AND BC_GRIDS

    !Local Variables
    !CLASS(POINT)	   :: Temp_Node ! Temp variable for subroutine compatibility. ! It's posible to use polymorphic variables
    INTEGER           :: i         ! Loop index

    DOUBLE PRECISION, DIMENSION(3) :: IndSpeed, xyz1, xyz2  ! Cumulative Induced Speed								!**
    INTEGER            :: N1, N2, k

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
    DOUBLE PRECISION, PARAMETER      :: PI=3.1415926535897933d0 ! Pi.
    DOUBLE PRECISION                 :: L               ! Length of the vortex segment.
    DOUBLE PRECISION                 :: delta    ! "Cutoff" radius parameter (VAN GARREL)
    END FUNCTION VtxSgmnt_IndSpeed
    end interface

     !$OMP PARALLEL SHARED (GRD, BodySet, nBodies)
     !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(i, N1, N2, k, IndSpeed, xyz1, xyz2)
    
    DO i=1, GRD%nPanels

        !Temp_Node%xyz = GRD%PanelSet(i)%CtrlPoint%xyz !Ctrl Point coords extraction

        !Temp_Node = GRD%PanelSet(i)%CtrlPoint
        call IndSpeed_atNode (BodySet, nBodies, GRD%PanelSet(i)%CtrlPoint) ! This subroutine in set for a polymorphic Variable.

        !Substracts the panel own influence
        IndSpeed = [0.d0, 0.d0, 0.d0]
        DO k=1, GRD%PanelSet(i)%nPanelSegments
            !	print *, IndSpeed
            N1 = GRD%VortexSegmentSet(abs(GRD%PanelSet(i)%PanelSegments(k)))%SegmentNodes(1)
            !	print *, 'N1: ',N1
            N2 = GRD%VortexSegmentSet(abs(GRD%PanelSet(i)%PanelSegments(k)))%SegmentNodes(2)
            !	print *, 'N2: ',N2
            xyz1 = GRD%NodalPointSet(N1)%xyz
            !	print *, 'xyz1: ',xyz1
            xyz2 = GRD%NodalPointSet(N2)%xyz
            !	print *, 'xyz2: ',xyz2
            !
            IndSpeed = IndSpeed + VtxSgmnt_IndSpeed (xyz1, xyz2, GRD%PanelSet(i)%CtrlPoint%xyz, GRD%VortexSegmentSet(abs(GRD%PanelSet(i)%PanelSegments(k)))%SgmntCirc)  !!! CHECKED: OK
            !	print *, IndSpeed
        END DO
        GRD%PanelSet(i)%CtrlPoint%InducedVelocity = GRD%PanelSet(i)%CtrlPoint%InducedVelocity - IndSpeed

        !call Wake_Influence_IndSpeed (BodySet(1)%Wakes(1), GRD%PanelSet(i)%CtrlPoint%xyz, WkeIndSpeed)  !**
        !GRD%PanelSet(i)%CtrlPoint%InducedVelocity = WkeIndSpeed ! Induced speed storage.				!**

    END DO
    
     !$OMP END DO
     !$OMP END PARALLEL

    END SUBROUTINE IndSpeed_atCtrlPoints
    
    
    !SUBROUTINE IndSpeed_atCtrlPoints_old (GRD, BodySet, nBodies)
    !
    !!Calculates the induced speed at every Control Point of the Grid (GRD).
    !!Considering all Wakes and Grids (BC_Grids as well) influence.
    !!Free stream speed (Vinf) is NOT considered.
    !! USES POLYMORPHIC VARIABLES SO AS TO PROCESS GRIDS AND BC_GRIDS.
    !
    !!Input Variables
    !type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    !INTEGER, INTENT(IN)                                 :: nBodies ! Number of bodies of the simulation.
    !
    !!Output Variables
    !CLASS(GRID), INTENT(INOUT)   :: GRD   !Provides the control points for induced speed calculation. POLYMORPHIC FOR GIRDS AND BC_GRIDS
    !
    !!Local Variables
    !!CLASS(POINT)	   :: Temp_Node ! Temp variable for subroutine compatibility. ! It's posible to use polymorphic variables
    !INTEGER           :: i         ! Loop index
    !
    !DOUBLE PRECISION, DIMENSION(3) :: IndSpeed, xyz1, xyz2  ! Cumulative Induced Speed								!**
    !INTEGER            :: N1, N2, k
    !
    !interface
    !FUNCTION VtxSgmnt_IndSpeed (NP1, NP2, CP, Gamma)
    !!Calculates the speed induced by a vortex segment (from NP1 to NP2) of intensity Gamma, in the control point CP.
    !
    !IMPLICIT NONE
    !
    !!Input Variables
    !DOUBLE PRECISION, DIMENSION (3)  :: NP1, NP2    ! Nodal Points of the Segment, first-last.
    !DOUBLE PRECISION, DIMENSION (3)  :: CP          ! Control Point to calculate induced speed in.
    !DOUBLE PRECISION                 :: Gamma       ! Circulation of the segment.
    !
    !!Output Variables
    !DOUBLE PRECISION, DIMENSION (3)  :: VtxSgmnt_IndSpeed   ! Induced speed in the control point.
    !
    !!Local Variables
    !DOUBLE PRECISION, DIMENSION (3)  :: R1, R2          ! Pointing vectors from CP to NPi.
    !DOUBLE PRECISION, PARAMETER      :: PI=3.1415926535 ! Pi.
    !DOUBLE PRECISION                 :: L               ! Length of the vortex segment.
    !DOUBLE PRECISION                 :: delta = 0.05    ! "Cutoff" radius parameter (VAN GARREL)
    !END FUNCTION VtxSgmnt_IndSpeed
    !end interface
    !
    !DO i=1, GRD%nPanels
    !
    !    !Temp_Node%xyz = GRD%PanelSet(i)%CtrlPoint%xyz !Ctrl Point coords extraction
    !
    !    !Temp_Node = GRD%PanelSet(i)%CtrlPoint
    !    call IndSpeed_atNode (BodySet, nBodies, GRD%PanelSet(i)%CtrlPoint) ! This subroutine in set for a polymorphic Variable.
    !
    !    !Substracts the panel own influence
    !    IndSpeed = [0.0, 0.0, 0.0]
    !    DO k=1, GRD%PanelSet(i)%nPanelSegments
    !        !	print *, IndSpeed
    !        N1 = GRD%VortexSegmentSet(abs(GRD%PanelSet(i)%PanelSegments(k)))%SegmentNodes(1)
    !        !	print *, 'N1: ',N1
    !        N2 = GRD%VortexSegmentSet(abs(GRD%PanelSet(i)%PanelSegments(k)))%SegmentNodes(2)
    !        !	print *, 'N2: ',N2
    !        xyz1 = GRD%NodalPointSet(N1)%xyz
    !        !	print *, 'xyz1: ',xyz1
    !        xyz2 = GRD%NodalPointSet(N2)%xyz
    !        !	print *, 'xyz2: ',xyz2
    !        !
    !        IndSpeed = IndSpeed + VtxSgmnt_IndSpeed (xyz1, xyz2, GRD%PanelSet(i)%CtrlPoint%xyz, GRD%VortexSegmentSet(abs(GRD%PanelSet(i)%PanelSegments(k)))%SgmntCirc)  !!! CHECKED: OK
    !        !	print *, IndSpeed
    !    END DO
    !    GRD%PanelSet(i)%CtrlPoint%InducedVelocity = GRD%PanelSet(i)%CtrlPoint%InducedVelocity - IndSpeed
    !
    !    !call Wake_Influence_IndSpeed (BodySet(1)%Wakes(1), GRD%PanelSet(i)%CtrlPoint%xyz, WkeIndSpeed)  !**
    !    !GRD%PanelSet(i)%CtrlPoint%InducedVelocity = WkeIndSpeed ! Induced speed storage.				!**
    !
    !END DO
    !
    !END SUBROUTINE IndSpeed_atCtrlPoints_old
    
    SUBROUTINE IndSpeed_atExtPoint (ExtPoint_Set, BodySet, nBodies)
    !Calculates the induced speed at every Extra Point of the Set.
    !Free stream speed (Vinf) is NOT considered.

    !Input Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    INTEGER, INTENT(IN)                                 :: nBodies ! Number of bodies of the simulation.

    !Output Variables
    type(Extra_Point), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: ExtPoint_Set

    !Local Variables
    INTEGER   :: i         ! Loop index
    INTEGER   :: NumberOfExtraPoints

    NumberOfExtraPoints = size(ExtPoint_Set)
    
     !$OMP PARALLEL SHARED (ExtPoint_Set, BodySet, nBodies)
     !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(i)
    
    DO i=1, NumberOfExtraPoints

        call IndSpeed_atNode (BodySet, nBodies, ExtPoint_Set(i)) ! This subroutine in set for a polymorphic Variable.

    END DO
    
     !$OMP END DO 
     !$OMP END PARALLEL
    
    
    END SUBROUTINE IndSpeed_atExtPoint

    !SUBROUTINE IndSpeed_atExtPoint_old (ExtPoint_Set, BodySet, nBodies)
    !!Calculates the induced speed at every Extra Point of the Set.
    !!Free stream speed (Vinf) is NOT considered.
    !
    !!Input Variables
    !type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    !INTEGER, INTENT(IN)                                 :: nBodies ! Number of bodies of the simulation.
    !
    !!Output Variables
    !type(Extra_Point), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)   :: ExtPoint_Set
    !
    !!Local Variables
    !INTEGER   :: i         ! Loop index
    !INTEGER   :: NumberOfExtraPoints
    !
    !NumberOfExtraPoints = size(ExtPoint_Set)
    !
    !DO i=1, NumberOfExtraPoints
    !
    !    call IndSpeed_atNode (BodySet, nBodies, ExtPoint_Set(i)) ! This subroutine in set for a polymorphic Variable.
    !
    !
    !END DO
    !
    !END SUBROUTINE IndSpeed_atExtPoint_old

    SUBROUTINE IndSpeed_atNode (BodySet, nBodies, Node)
    ! Calculates the induced speed at a specific Nodal Point (Node) produced by grids and wakes influence.
    ! USES POLYMORPHIC VARIABLES SO AS TO PROCESS NODAL_POINTS AND CONTROL_POINTS.

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
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.
    INTEGER, INTENT(IN)                                 :: nBodies ! Number of bodies of the simulation.

    !Output Variables
    CLASS(POINT), INTENT(INOUT) :: Node    !Node where the induced speed is calculated. POLYMORPHIC FOR NODAL_POINTS \ CONTROL_POINTS \ EXTRA_POINTS

    !type(NODAL_POINT), INTENT(INOUT) :: Node    !Node where the induced speed is calculated

    !LocalVariables
    INTEGER             :: i,j,k  !Loop Indexes
    DOUBLE PRECISION, DIMENSION (3) :: IndSpeed, WkeIndSpeed  ! Cumulative induced speed
    INTEGER             :: N1,N2   ! Temp nodes id (Label=Position)
    INTEGER             :: L1, L2  ! Temp nodes lables
    DOUBLE PRECISION, DIMENSION (3) :: xyz1, xyz2, NodeXyz ! Temp nodes coordinates
    INTEGER             :: NodesStride !Size of the block of nodes of the wake (=InitNodes)
    INTEGER             :: CNT

    IndSpeed =  [0.0,0.0,0.0]
    !WkeIndSpeed =  [0.0,0.0,0.0]

    !NodesStride = WKE%NodesStride

    ! Grids Influence
    DO i=1,nBodies
        DO j=1,BodySet(i)%nGrids
            Grdj: associate (GRD => BodySet(i)%Grids(j))
                GrdVS: associate (nGrdVS => GRD%nVortexSegments)

                    DO k=1,nGrdVS
                        N1 = GRD%VortexSegmentSet(k)%SegmentNodes(1)
                        N2 = GRD%VortexSegmentSet(k)%SegmentNodes(2)
                        xyz1 = GRD%NodalPointSet(N1)%xyz
                        xyz2 = GRD%NodalPointSet(N2)%xyz

                        IndSpeed = IndSpeed + VtxSgmnt_IndSpeed (xyz1, xyz2, Node%xyz, GRD%VortexSegmentSet(k)%SgmntCirc)  !!! CHECKED: OK

                    END DO

                end associate GrdVS
            end associate Grdj

        END DO
    END DO

    ! BC_Grids Influence
    DO i=1,nBodies
        DO j=1,BodySet(i)%nBCGrids
            BCGrdj: associate (BCGRD => BodySet(i)%BCGrids(j))
                GrdVS: associate (nGrdVS => BCGRD%nVortexSegments)

                    DO k=1,nGrdVS
                        N1 = BCGRD%VortexSegmentSet(k)%SegmentNodes(1)
                        N2 = BCGRD%VortexSegmentSet(k)%SegmentNodes(2)
                        xyz1 = BCGRD%NodalPointSet(N1)%xyz
                        xyz2 = BCGRD%NodalPointSet(N2)%xyz

                        IndSpeed = IndSpeed + VtxSgmnt_IndSpeed (xyz1, xyz2, Node%xyz, BCGRD%VortexSegmentSet(k)%SgmntCirc)  !!! CHECKED: OK

                    END DO

                end associate GrdVS
            end associate BCGrdj

        END DO
    END DO
    
        ! print*, 'GridInducedSp >>  ',IndSpeed
    ! print*, 'WakeInducedSp >>  ',WkeIndSpeed
select type (Node) ! Type selector: consider the case of processing a BC_GRID in which the transporation velocity must be considered in the RHS.
type is (Nodal_Point)

    ! Wakess Influence
    DO i=1,nBodies
        DO j=1,BodySet(i)%nWakes
            Wkej: associate (WKE => BodySet(i)%Wakes(j))

                call Wake_Influence_IndSpeed (WKE, Node%xyz, WkeIndSpeed)

                IndSpeed = IndSpeed + WkeIndSpeed

            end associate Wkej

        END DO
    END DO
        Node%InducedSpeed = IndSpeed

type is (Control_Point)
		
		Node%InducedVelocity = IndSpeed + Node%WKEInducedVelocity  ! The wakes' induced speed was already calculated for RHS

type is (Extra_Point)
		
			! Wakess Influence
		DO i=1,nBodies
		    DO j=1,BodySet(i)%nWakes
		        Wkej: associate (WKE => BodySet(i)%Wakes(j))

		            call Wake_Influence_IndSpeed (WKE, Node%xyz, WkeIndSpeed)
		            IndSpeed = IndSpeed + WkeIndSpeed

		        end associate Wkej

		    END DO
		END DO
	
        Node%Velocity = IndSpeed
    end select

    END SUBROUTINE IndSpeed_atNode

SUBROUTINE DGDt_calc (Pnl, TStp, DGDt)
    !Calculates the DeltaCp unsteady component DG(t)/Dt for a panel. Approximates by a first-order finite difference.

    !Input Variables
    type(PANEL), INTENT(IN) :: Pnl
    DOUBLE PRECISION, INTENT(IN)        :: TStp

    !Output Variables
    DOUBLE PRECISION, INTENT(OUT)       :: DGDt

    DGDt = (Pnl%RingVorticity - Pnl%RngVtctyPrvStp) / TStp


END SUBROUTINE DGDt_calc

    SUBROUTINE DeltaCp_calc (BodySet, nBodies, TStp, Vinf, StepCnt)
    !Calculates the DeltaCp for each panel of a Grid or Bc_grid. Stores it in the Control Point Variable.
    !Uses Grids/BC_Grids and Wakes from each body,and the time step variable.

    !Input Variables
    DOUBLE PRECISION, INTENT(IN)                :: TStp  ! Time step
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: Vinf  ! Free stream velocity
    INTEGER, INTENT(IN)             :: nBodies ! Number of bodies of the simulation
    INTEGER, INTENT(IN)             :: StepCnt ! Step Number

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation. Provides Grids and Wakes data.

    !Local Variables
    DOUBLE PRECISION, DIMENSION(3)  :: DeltaV ! Change in tangencial velocity across the panel
    DOUBLE PRECISION, DIMENSION(3)  :: OwnTspVel ! Aux variable for extracting OwnVelocity and adding TranspVelocity when BC_Grid is processed.
    DOUBLE PRECISION                :: DGDt ! Unsteady change in panel ring vorticity term.
    DOUBLE PRECISION                :: DCp !Delta Cp temp variable
    INTEGER             :: i,j,k !Loop indexes


    ! Grids

    DO i=1,nBodies
        DO j=1,BodySet(i)%nGrids
            Grdj: associate (GRD => BodySet(i)%Grids(j))

                call Area_calc (GRD%PanelSet, GRD%VortexSegmentSet, GRD%NodalPointSet) ! Calculates the area of each panel.
                !call IndSpeed_atCtrlPoints (GRD, BodySet, nBodies) ! Calculates induced speed at control points.

                DO k=1, GRD%nPanels
                    call Panel_DeltaV (GRD%PanelSet(k), GRD%VortexSegmentSet, GRD%NodalPointSet, DeltaV)
                    call DGDt_calc (GRD%PanelSet(k), TStp, DGDt)
                    !print *, 'DGDt: ', DGDt

                    !if (StepCnt == 1) then
                    !	DGDt = 0.0
                    !end if

                    OwnTspVel= GRD%PanelSet(k)%CtrlPoint%OwnVelocity

                    DCp = 2.d0*(dot_product((GRD%PanelSet(k)%CtrlPoint%InducedVelocity + Vinf),DeltaV) - DGDt + dot_product(DeltaV,-OwnTspVel))!**********-*-*-*-*--*-*-*-*-*-*

                    !print*, 'indVel: ', GRD%PanelSet(k)%CtrlPoint%InducedVelocity
                    !print*, 'deltaV: ', DeltaV
                    !pause

                    !DCp = DCp * (sign(1.0,GRD%PanelSet(k)%RingVorticity))

                    GRD%PanelSet(k)%CtrlPoint%Delta_Cp = DCp / norm2(Vinf)**2.d0 !!!! Adimensionalization
                    !GRD%PanelSet(k)%CtrlPoint%Delta_Cp = DCp / (77.8)**2 !!!! Adimensionalization
                    !GRD%PanelSet(k)%CtrlPoint%Delta_Cp = DCp / (36156.d0)**2 !!!! Adimensionalization
                END DO
                !pause
            end associate Grdj
        END DO
    END DO

    ! BC_Grids  NOT DEFINED :::: Sets equal to 1.0
    DO i=1,nBodies
        DO j=1,BodySet(i)%nBCGrids
            BcGrdj: associate (BCGRD => BodySet(i)%BCGrids(j))

                !call Area_calc (BCGRD%PanelSet, BCGRD%VortexSegmentSet, BCGRD%NodalPointSet) ! Calculates the area of each panel.
                !!call IndSpeed_atCtrlPoints (BCGRD, BodySet, nBodies) ! Calculates induced speed at control points.

                DO k=1, BCGRD%nPanels

                    !call Panel_DeltaV (BCGRD%PanelSet(k), BCGRD%VortexSegmentSet, BCGRD%NodalPointSet, DeltaV)
                    !call DGDt_calc (BCGRD%PanelSet(k), TStp, DGDt)
                    !
                    !if (BCGRD%TrnspCondition==1) then  ! Case for Transp velocity constant (non-zero) over the grid.
                    !    OwnTspVel= BCGRD%PanelSet(k)%CtrlPoint%OwnVelocity + BCGRD%TranspirationVel * BCGRD%PanelSet(k)%CtrlPoint%NormalVersor
                    !else if (BCGRD%TrnspCondition==2) then ! Case for Transp velocity variable for each panel.
                    !    OwnTspVel= BCGRD%PanelSet(k)%CtrlPoint%OwnVelocity + BCGRD%PanelSet(k)%CtrlPoint%TrnspVelocity * BCGRD%PanelSet(k)%CtrlPoint%NormalVersor
                    !else ! Case BC_GRID with no penetration condition
                    !    OwnTspVel= BCGRD%PanelSet(k)%CtrlPoint%OwnVelocity
                    !end if
                    !
                    !DCp = 2.0*(dot_product((BCGRD%PanelSet(k)%CtrlPoint%InducedVelocity + Vinf),DeltaV) - (DGDt + dot_product(DeltaV,OwnTspVel)))
                    !!DCp = DCp * (sign(1.0,BCGRD%PanelSet(k)%RingVorticity))
                    !
                    !BCGRD%PanelSet(k)%CtrlPoint%Delta_Cp = DCp / norm2(Vinf)**2 !!!! Adimensionalization
                    !
                    BCGRD%PanelSet(k)%CtrlPoint%Delta_Cp = 1.d0

                END DO
            end associate BcGrdj
        END DO
    END DO

    END SUBROUTINE DeltaCp_calc

SUBROUTINE BC_Segment_Condition (BC_GRD)
    !Determines the condition of a BC_grid Vortex Segment respecting the "edge condition", i.e. "Inner Segment", "Non-Convective Edge".
    !For DeltaCP calculation (if needed)

    !Input Variables
    !None

    !Output Variables
    type(BC_GRID), INTENT(INOUT)         :: BC_GRD

    !Local Variables
    INTEGER :: i, j !Loop indexes
    INTEGER :: GrdVSLabel, ClnVSLabel ! Labels extraction variables

    DO i=1,BC_GRD%nVortexSegments
        BC_GRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 1.d0
        BC_GRD%VortexSegmentSet(i)%Sgmnt_Cond(2) = 0.d0

        !Edge Identification
        if (BC_GRD%VortexSegmentSet(i)%nSharingPanels == 1) then
            BC_GRD%VortexSegmentSet(i)%Sgmnt_Cond(1) = 2.d0
        end if
    END DO

END SUBROUTINE BC_Segment_Condition

SUBROUTINE ForceCoeficient_calc (GRD, V_inf, CL, CD, CM, Ref_Pt, Ref_Vrs)
    !Calculates the CL and CD for the grid at every step.

    !Input Variables
    type(GRID), INTENT(IN)         :: GRD
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: V_inf
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Ref_Pt  ! Reference axis point
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Ref_Vrs ! Reference axis versor

    !Output Variables
    DOUBLE PRECISION, INTENT(OUT) :: CL
    DOUBLE PRECISION, INTENT(OUT) :: CD
    DOUBLE PRECISION, INTENT(OUT) :: CM

    !Local Variables
    INTEGER            :: i                              ! Loop Index
    DOUBLE PRECISION, DIMENSION(3) :: D_Cp, Mmnt                     ! Delta CP cummulative variable, Moment cummulative variable
    DOUBLE PRECISION, DIMENSION(3) :: CFv                            ! Aux Force Coeficient vectorial variable
    DOUBLE PRECISION, DIMENSION(3) :: u_vect, s_vect, r_vect, F_vect ! Aux vectors
    DOUBLE PRECISION               :: Area, Area_1                   ! Area cummulative variable, 1/Area

    D_Cp = [0.0,0.0,0.0]
    Area = 0.0
    Mmnt = [0.0,0.0,0.0]

    DO i=1,GRD%nPanels
        !Force Coeficient
        D_Cp = D_Cp + GRD%PanelSet(i)%CtrlPoint%Delta_Cp * GRD%PanelSet(i)%Area * GRD%PanelSet(i)%CtrlPoint%NormalVersor
        Area = Area + GRD%PanelSet(i)%Area

        !!Moment Coeficient
        !u_vect = GRD%PanelSet(i)%CtrlPoint%xyz - Ref_Pt
        !s_vect = dot_product(u_vect,Ref_Vrs) * Ref_Vrs
        !r_vect = u_vect - s_vect

        !print*,'u = ', u_vect
        !print*,'s = ', s_vect
        !print*,'r = ', r_vect

        F_vect = GRD%PanelSet(i)%CtrlPoint%Delta_Cp * GRD%PanelSet(i)%Area * GRD%PanelSet(i)%CtrlPoint%NormalVersor
        !Mmnt = Mmnt + cross_product(r_vect,F_vect)

    END DO

    print*,'AREA = ', Area
    print*,'D_Cp = ', D_Cp
    !print*,'Mmnt = ', Mmnt
    !pause

    Area_1 = 1.d0/Area

    CFv = D_Cp * Area_1
    print*,'CFv = ', CFv

    CD = dot_product(CFv,V_inf*(1.0/norm2(V_inf)))

    !CL= norm2(CFv - CD*V_inf*(1/norm2(V_inf)))

    CL = sqrt(norm2(CFv)**2.0 - CD**2.0) * sign(1.0,CFv(3)) !Uses (3) for z-component, change otherwise

    CM = norm2(CFv)
    !CM = dot_product(Mmnt,Ref_Vrs)*Area_1  ! Moment Coefficient*Chord   Must be divided by the chord

    !CL = norm2(CFv) !!! ONLY FOR ROTATING BLADES
    !CD = CFv(2)  ! Val check

    !print*,'CD = ', CD
    !print*,'CL = ', CL
    !print*,'CM = ', CM

    !pause

    !pause
END SUBROUTINE ForceCoeficient_calc


END MODULE Aero_Loads