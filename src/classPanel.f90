    MODULE classPanel

    ! "Panel" abstract data type definition

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE classControl_Point
    USE classVortex_Segment
    USE classNodal_Point, only: Nodal_Point
    USE ShFunct_DShFunct
    USE P_S_N

    IMPLICIT NONE
    PUBLIC

    TYPE, Public :: Panel
        !Attributes:
        INTEGER                             :: Grid ! Owner
        INTEGER                             :: Label ! Global Id number of the panel
        INTEGER                             :: Panel_Type     ! Sets the panel configuration according to internal code.
        INTEGER                             :: nPanelNodes    ! Number of nodes of the panel (3:8)
        INTEGER, ALLOCATABLE, DIMENSION (:) :: PanelNodes     ! Nodes of the panel (3:8)
        INTEGER                             :: nPanelSegments ! Number of segments of the panel, matches the mumber of nodes.
        INTEGER, ALLOCATABLE, DIMENSION (:) :: PanelSegments  ! vortex segments of the panel
        DOUBLE PRECISION                    :: RingVorticity  ! vortex ring vorticity of the panel
        DOUBLE PRECISION                    :: RngVtctyPrvStp ! vortex ring vorticity of the panel in the previous time step (For deltaCp calculation)
        !DOUBLE PRECISION, DIMENSION(3)                  :: DeltaV        ! Jump in tangential speed across the panel. (For deltaCp calculation)
        type(CONTROL_POINT)                 :: CtrlPoint
        DOUBLE PRECISION                    :: Area           ! Area of the panel

    CONTAINS ! methods' names

    PROCEDURE, Public, NoPass :: PanelSet_Allocation
    PROCEDURE, Public, NoPass :: PanelSet_Generation
    PROCEDURE, Public, NoPass :: PanelSet_Printing
    PROCEDURE, Public, NoPass :: Xa_Calculation
    PROCEDURE, Public, NoPass :: Sgmnt_Coef
    PROCEDURE, Public, NoPass :: Panel_DeltaV
    PROCEDURE, Public, NoPass :: Area_calc

    END TYPE Panel

    CONTAINS ! ==================================================================

    SUBROUTINE PanelSet_Allocation (NumberOfPanels, Pnl_Set)  ! Set allocation (Panel allocatable array)

    INTEGER, INTENT(IN)                                    :: NumberOfPanels
    type(PANEL), ALLOCATABLE, DIMENSION (:), INTENT(INOUT) :: Pnl_Set ! Array set variable for allocation

    ALLOCATE (Pnl_Set (NumberOfPanels))

    END SUBROUTINE PanelSet_Allocation

    SUBROUTINE  PanelSet_Generation (GridLabel, IENArray, IESArray, NodesArray, PanelType, PanelSet, V_inf)
    ! Generates a Grid's PanelSet from the input data arrays, works with P_S_N Module.

    !Input Variables
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(IN) :: IENArray
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(IN) :: IESArray
    DOUBLE PRECISION, DIMENSION (:,:), INTENT(IN)                 :: NodesArray
    INTEGER, INTENT(IN)                               :: GridLabel
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)                    :: V_inf !Free stream speed.
    INTEGER, DIMENSION (:), ALLOCATABLE, INTENT(IN)   :: PanelType

    !Output Variables
    type(PANEL), DIMENSION (:), INTENT(INOUT)  :: PanelSet

    !Local Variables
    INTEGER               :: i, j, NumberOfPanels, nP_NS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Xa
    DOUBLE PRECISION, DIMENSION (8)   :: PShapeFunct
    DOUBLE PRECISION, DIMENSION (8,2) :: Derived_PSF
    DOUBLE PRECISION, DIMENSION (3)   :: NrmlVersor, XIv, CP_xyz
    !INTEGER, DIMENSION (:), ALLOCATABLE  :: PanelType

    DOUBLE PRECISION,DIMENSION(3) :: OwnV = [0.d0,0.d0,0.d0] !|--------------------->>>> FOR TRIAL ONLY


    NumberOfPanels=size(PanelSet)

    do i=1,NumberOfPanels
        !LABEL
        PanelSet(i)%Grid = GridLabel
        PanelSet(i)%Label = i

        !PANEL TYPE
        PanelSet(i)%Panel_Type=PanelType(i)  !Defines where is each node in the panel, according to internal code

        select case (PanelType(i))  !Determines number of nodes (and of segments) of the panel from panel type
        case (1)
            nP_NS = 4
        case (3)
            nP_NS = 3
        case (51, 52, 53, 54)
            nP_NS = 5
        case (61, 62, 63, 64)
            nP_NS = 6
        case (71, 72, 73, 74)
            nP_NS = 7
        case (8)
            nP_NS = 8
        end select

        !NODES
        PanelSet(i)%nPanelNodes = nP_NS
        ALLOCATE (PanelSet(i)%PanelNodes (PanelSet(i)%nPanelNodes))
        do j=1,4   !--------------------------------------------MODIFIED FOR TRIAL Original: do j=1,8
            if (IENArray(j,i) /= 0) then   !Excludes IEN elements with zero value => PanelType must be kept
                PanelSet(i)%PanelNodes(j) = IENArray(j,i)  ! CHECK THE INDEX!!!
            end if
        end do

        !SEGMENTS
        PanelSet(i)%nPanelSegments= nP_NS
        ALLOCATE (PanelSet(i)%PanelSegments (PanelSet(i)%nPanelSegments))
        do j=1,8
            if (IENArray(j,i) /= 0) then   !Excludes IES elements with zero value => PanelType must be kept
                PanelSet(i)%PanelSegments= IESArray(:,i)
            end if
        end do

        !CONTROL POINT
        XIv=[0,0,0]
        call Xa_Calculation (PanelSet(i)%PanelNodes, NodesArray, PanelSet(i)%Panel_Type, Xa)
        call Shape_Functions (PanelSet(i)%Panel_Type, PShapeFunct, XIv)
        !print *, 'Ptype', PanelSet(i)%Panel_Type
        call Derived_Shape_Functions (PanelSet(i)%Panel_Type, Derived_PSF)
        !print *, 'DPSF', Derived_PSF
        call CP_Coords (Xa, PShapeFunct, CP_xyz)
        call CP_NormalVersor (Xa, Derived_PSF, NrmlVersor)
        !print *, 'xyz', CP_xyz
        PanelSet(i)%CtrlPoint%xyz=CP_xyz
        PanelSet(i)%CtrlPoint%xyz_loc=CP_xyz
        !print *, 'nor', NrmlVersor
        PanelSet(i)%CtrlPoint%NormalVersor=NrmlVersor
        !----------------------------------------------------------+
        call CP_Velocity (PanelSet(i)%CtrlPoint, V_inf, OwnV)      !| -------> for trial only
        !----------------------------------------------------------+

        !VORTICITY INITIALIZATION
        PanelSet(i)%RingVorticity = 0.d0
        PanelSet(i)%RngVtctyPrvStp = 0.d0
    end do

    END SUBROUTINE PanelSet_Generation

    SUBROUTINE PanelSet_Printing (Pnl_Set, FileUnit) ! Prints a Grid's PanelSet
    !File unit specification is needed

    !Input Variables
    INTEGER, INTENT(IN)                     :: FileUnit  !File id number for writing
    type(PANEL), DIMENSION (:), INTENT(IN)  ::  Pnl_Set ! Array set variable for Printing

    !Output Variables
    !None

    !Local Variables
    INTEGER             :: NumberOfPanels,i

    NumberOfPanels=size(Pnl_Set)

    !write (FileUnit,'(A16,/)') 'Panel_Set Data: '

    !write (FileUnit,'(A82,/)') 'Grid:  | Label: |	  	  Nodes:	       |						Segments:		                | Area: |'

    !do i=1,NumberOfPanels
    !    write (FileUnit, 60) Pnl_Set(i)%Grid,Pnl_Set(i)%Label,Pnl_Set(i)%PanelNodes,Pnl_Set(i)%PanelSegments,Pnl_Set(i)%Area
    !end do
    !Ctrl Point data
    !write (FileUnit,'(A160,/)') 'Label: |      Coords:         |          Own:            |     Induced:             |         Wind:             |          Normal:         |   CtrlPt_DCp   |'
    !do i=1,NumberOfPanels
    !     write (FileUnit, 61) Pnl_Set(i)%Label,Pnl_Set(i)%CtrlPoint%xyz,Pnl_Set(i)%CtrlPoint%OwnVelocity,Pnl_Set(i)%CtrlPoint%InducedVelocity,Pnl_Set(i)%CtrlPoint%WKEInducedVelocity,Pnl_Set(i)%CtrlPoint%WindVelocity,Pnl_Set(i)%CtrlPoint%NormalVersor,'   ',Pnl_Set(i)%CtrlPoint%Delta_Cp
    !end do
    
    do i=1,NumberOfPanels
        write (FileUnit, 61),Pnl_Set(i)%CtrlPoint%Delta_Cp,'   '
	end do

	
60  FORMAT (2I7,4I7,8I7,F9.3)    !Space for more nodes or segments must be added.
61  FORMAT (F16.12,A3,$)

	END SUBROUTINE PanelSet_Printing
	
!	SUBROUTINE PanelSet_Printing (Pnl_Set, FileUnit) ! Prints a Grid's PanelSet
!    !File unit specification is needed
!
!    !Input Variables
!    INTEGER, INTENT(IN)                     :: FileUnit  !File id number for writing
!    type(PANEL), DIMENSION (:), INTENT(IN)  ::  Pnl_Set ! Array set variable for Printing
!
!    !Output Variables
!    !None
!
!    !Local Variables
!    INTEGER             :: NumberOfPanels,i
!
!    NumberOfPanels=size(Pnl_Set)
!
!    write (FileUnit,'(A16,/)') 'Panel_Set Data: '
!
!    write (FileUnit,'(A82,/)') 'Grid:  | Label: |	  	  Nodes:	       |						Segments:		                | Area: |'
!
!    !do i=1,NumberOfPanels
!    !    write (FileUnit, 60) Pnl_Set(i)%Grid,Pnl_Set(i)%Label,Pnl_Set(i)%PanelNodes,Pnl_Set(i)%PanelSegments,Pnl_Set(i)%Area
!    !end do
!    !Ctrl Point data
!    write (FileUnit,'(A160,/)') 'Label: |      Coords:         |          Own:            |     Induced:             |         Wind:             |          Normal:         |   CtrlPt_DCp   |'
!    do i=1,NumberOfPanels
!        write (FileUnit, 61) Pnl_Set(i)%Label,Pnl_Set(i)%CtrlPoint%xyz,Pnl_Set(i)%CtrlPoint%OwnVelocity,&
!            Pnl_Set(i)%CtrlPoint%InducedVelocity,Pnl_Set(i)%CtrlPoint%WindVelocity,Pnl_Set(i)%CtrlPoint%NormalVersor,'   ',Pnl_Set(i)%CtrlPoint%Delta_Cp
!    end do
!
!60  FORMAT (2I7,4I7,8I7,F9.3)    !Space for more nodes or segments must be added.
!61  FORMAT (I3,3F9.3,3F9.3,3F9.3,3F9.3,3F9.3,A3,F16.12)
!
!    END SUBROUTINE PanelSet_Printing

    SUBROUTINE Xa_Calculation (PanelNodes, NodesArray, PanelType, Xa)
    ! Calculates a matrix with the coordinates of each node of a panel as columns, according to the Panel_Type configuration.

    !Input Variables
    INTEGER, DIMENSION (:), INTENT(IN) :: PanelNodes ! Panel Nodes ID
    DOUBLE PRECISION, DIMENSION (:,:), INTENT(IN)  :: NodesArray ! Array of Nodes' coordinates
    INTEGER, INTENT (IN)               :: PanelType  ! Sets the configuration of the Xa matrix.

    !Output Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT) :: Xa     ! Output matrix

    !Local Variables
    INTEGER :: nNodes ! Number of nodes of the panel - ID of each node
    INTEGER :: i ! Loop index

    nNodes = size(PanelNodes)

    ALLOCATE (Xa(3,8))
    Xa = 0 ! Initialize with zeros.

    if (PanelType == 3) then
        do i=1,3    !CASE 3.
            Xa(:,i) =  NodesArray(:,PanelNodes(i))
        end do
    else
        do i=1,4   !CASE 1. Applies for all except for CASE 3.
            Xa(:,i) =  NodesArray(:,PanelNodes(i))
        end do
    end if

    select case (PanelType)

    case (51)
        Xa(:,5) =  NodesArray(:,PanelNodes(5))
    case (52)
        Xa(:,6) =  NodesArray(:,PanelNodes(6))
    case (53)
        Xa(:,7) =  NodesArray(:,PanelNodes(7))
    case (54)
        Xa(:,8) =  NodesArray(:,PanelNodes(8))
    case (61)
        Xa(:,5) =  NodesArray(:,PanelNodes(5))
        Xa(:,6) =  NodesArray(:,PanelNodes(6))
    case (62)
        Xa(:,6) =  NodesArray(:,PanelNodes(6))
        Xa(:,7) =  NodesArray(:,PanelNodes(7))
    case (63)
        Xa(:,7) =  NodesArray(:,PanelNodes(7))
        Xa(:,8) =  NodesArray(:,PanelNodes(8))
    case (64)
        Xa(:,5) =  NodesArray(:,PanelNodes(5))
        Xa(:,8) =  NodesArray(:,PanelNodes(8))
    case (71)
        Xa(:,5) =  NodesArray(:,PanelNodes(5))
        Xa(:,6) =  NodesArray(:,PanelNodes(6))
        Xa(:,7) =  NodesArray(:,PanelNodes(7))
    case (72)
        Xa(:,6) =  NodesArray(:,PanelNodes(6))
        Xa(:,7) =  NodesArray(:,PanelNodes(7))
        Xa(:,8) =  NodesArray(:,PanelNodes(8))
    case (73)
        Xa(:,5) =  NodesArray(:,PanelNodes(5))
        Xa(:,7) =  NodesArray(:,PanelNodes(7))
        Xa(:,8) =  NodesArray(:,PanelNodes(8))
    case (74)
        Xa(:,5) =  NodesArray(:,PanelNodes(5))
        Xa(:,6) =  NodesArray(:,PanelNodes(6))
        Xa(:,8) =  NodesArray(:,PanelNodes(8))
    case (8)
        do i=5,8
            Xa(:,i) =  NodesArray(:,PanelNodes(i))
        end do

    end select

    END SUBROUTINE Xa_Calculation

    SUBROUTINE Sgmnt_Coef (Pnl, VtxSgmntSet, NodesSet)
    ! Generates a structure with the coeficients (+1,-1) for each vortex segment in reference to the
    ! ring vorticity determinated by the direction of the normal versor of each panel, according to the dextrorotatory rule.
    ! The structure is applied on the segment list of each panel adding a negative sign to the segment label if it is opposite
    ! to the ring vorticity direction defined. For example, the list PanelSegments=[14,-214,76,-99] indicates that the directions
    ! of segments #214 and #99 are contrary to the ring vorticity.


    interface
    FUNCTION cross_product(U,V)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
    DOUBLE PRECISION, DIMENSION(3)              :: cross_product
    DOUBLE PRECISION                            :: cross1, cross2, cross3

    END FUNCTION cross_product
    end interface

    ! Input Variables
    type(NODAL_POINT), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: NodesSet
    type(VORTEX_SEGMENT), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: VtxSgmntSet

    ! Output Variables
    type(PANEL), INTENT(INOUT) :: Pnl

    ! Local Variables
    INTEGER  :: i, nNodes, nSgmnts, N1, N2, SgmntLbl, alpha
    DOUBLE PRECISION, DIMENSION(3) :: N1coor, N2coor, Gamma, G2, n

    nNodes=size(NodesSet)
    nSgmnts=size(VtxSgmntSet)

    do i=1,Pnl%nPanelSegments
        SgmntLbl=Pnl%PanelSegments(i)     !Segment label
        N1=VtxSgmntSet(SgmntLbl)%SegmentNodes(1)    !Label of First node of the segment
        N2=VtxSgmntSet(SgmntLbl)%SegmentNodes(2)    !Label of Second node of the segment
        N1coor=NodesSet(N1)%xyz         !Coords of First node of the segment
        N2coor=NodesSet(N2)%xyz         !Coords of Second node of the segment
        Gamma=N2coor-N1coor         !Segment as vector
        G2=N2coor-Pnl%CtrlPoint%xyz       !Vector from panel CtrlPoint to Second node of the segment
        n=Pnl%CtrlPoint%NormalVersor      !Panel normal versor

        alpha = sign(1.0,dot_product(G2,cross_product(Gamma,n)))     !Segment Coeficient

        Pnl%PanelSegments(i)=Pnl%PanelSegments(i)*alpha

    end do

    END Subroutine Sgmnt_Coef

    SUBROUTINE Panel_DeltaV (Pnl, VtxSgmntSet, NodesSet, DeltaV)
    ! Calculates the change in the tangential speed across the panel, taking as inputs the corresponding VtxSgmnts and Nodes Sets.
    ! Applies for a grid structure.
    ! IF SEVEDOUBLE PRECISION GRIDS ARE CONSIDERED THE SUBROUTINE MUST BE MODIFIED SO AS TO INCLUDE THEM.

    interface
    FUNCTION cross_product(U,V)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
    DOUBLE PRECISION, DIMENSION(3)              :: cross_product
    DOUBLE PRECISION                            :: cross1, cross2, cross3

    END FUNCTION cross_product
    end interface

    !Input Variables
    type(Vortex_Segment), DIMENSION (:), INTENT(IN)  ::  VtxSgmntSet ! Corresponding Vtx Segments set
    type(Nodal_Point), DIMENSION (:), INTENT(IN)     ::  NodesSet    ! Corresponding Nodes set
    type(Panel), INTENT(IN)  ::  Pnl !For area value

    !Output Variables
    DOUBLE PRECISION, DIMENSION(3) :: DeltaV  ! Change in the tangential speed across the panel

    !Local Variables
    DOUBLE PRECISION               :: Area  ! Temp Variable. Area of the panel
    DOUBLE PRECISION, DIMENSION(3) :: Gamma ! Gamma vector definition
    DOUBLE PRECISION               :: Gc    ! Current Segments Circ extraction variable
    DOUBLE PRECISION               :: Gp    ! Previous Segments Circ extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: L   ! Segment extraction variable
    !DOUBLE PRECISION               :: Ln
    INTEGER            :: Nid1, Nid2 ! Nodes id extraction variable
    DOUBLE PRECISION               :: a, b  ! Sgmnt Condition coeficients extraction variable
	DOUBLE PRECISION               :: alpha ! Segment/panel area ratio: Area/SharedArea^2
    INTEGER            :: i     !Loop index
    !DOUBLE PRECISION               :: beta ! Convection Coeficient, for different sized panels.

    !if (Pnl%Panel_Type == 1) then ! SET ALGORITHM FOR DIFFERENT PANEL_TYPES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !Gamma vector
    !   Gamma = 0.0
    !
    !   DO i=1,Pnl%nPanelSegments
    !
    !      Nid1 = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SegmentNodes(1) ! Initial Node ID
    !      Nid2 = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SegmentNodes(2) ! Final Node ID
    !
    !      L = NodesSet(Nid2)%xyz - NodesSet(Nid1)%xyz ! Segment coords
    !
    !      a = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%Sgmnt_Cond(1) ! Segment condition coeficient "a"
    !      b = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%Sgmnt_Cond(2) ! Segment condition coeficient "b"
    !
    !!Check add
    !!select case (VtxSgmntSet(abs(Pnl%PanelSegments(i)))%Label)
    !!	case (13,10,7,4,238,236,234,231)
    !!
    !!		a = 1.0
    !!          print *, 'a : ', VtxSgmntSet(abs(Pnl%PanelSegments(i)))%Label, a
    !!end select
    !
    !      Gc = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SgmntCirc ! Actual Sgmnt Circulation
    !      Gp = sign(1,Pnl%PanelSegments(i))*Pnl%RngVtctyPrvStp    ! Previous Sgmnt Circulation
    !
    !      Gamma = Gamma + 0.5*(a)*L*(Gc-Gp*b)
    !
    !   END DO
    !DeltaV = -cross_product(Pnl%CtrlPoint%NormalVersor,Gamma)*(1/Pnl%Area)
    !end if
    !
    ! if (Pnl%Panel_Type == 3) then ! SET ALGORITHM FOR DIFFERENT PANEL_TYPES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !Gamma vector  ----
    Gamma = 0.0

    DO i=1,Pnl%nPanelSegments

        Nid1 = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SegmentNodes(1) ! Initial Node ID
        Nid2 = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SegmentNodes(2) ! Final Node ID

        L = NodesSet(Nid2)%xyz - NodesSet(Nid1)%xyz ! Segment coords

        !a = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%Sgmnt_Cond(1) ! Segment condition coeficient "a" ******
        b = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%Sgmnt_Cond(2) ! Segment condition coeficient "b"

        !print *, 'a,b: ',a,b

        Gc = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SgmntCirc ! Actual Sgmnt Circulation

        Gp = sign(1,Pnl%PanelSegments(i))*Pnl%RngVtctyPrvStp   ! Previous Sgmnt Circulation
        !beta = VtxSgmntSet(abs(Pnl%PanelSegments(i)))%Length/sqrt(VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SharedArea)
        !print *, 'Gp : ',Gp
        !Gamma = Gamma + 0.5*a*L*(Gc-Gp*b)!*********
		
		!alpha = Pnl%Area/(VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SharedArea*VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SharedArea*(b + 1.d0))

        Gamma = Gamma + L*(Gc-Gp*b)/ (VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SharedArea*(b + 1.d0))
		!Gamma = Gamma + 2.d0 * L*(Gc-Gp*b) * alpha
        !Gamma = Gamma + L*(Gc)/ (VtxSgmntSet(abs(Pnl%PanelSegments(i)))%SharedArea*(b+1))!

    END DO

    DeltaV = -cross_product(Pnl%CtrlPoint%NormalVersor,Gamma)!******
    !DeltaV = cross_product(Pnl%CtrlPoint%NormalVersor,Gamma)*(-1/Pnl%Area)  -----

    !end if

    END SUBROUTINE Panel_DeltaV

    SUBROUTINE Area_calc (PnlSet, VtxSgmntSet, NodesSet)
    !Calculates the area of each panel in the grid.

    interface
    FUNCTION cross_product(U,V)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
    DOUBLE PRECISION, DIMENSION(3)              :: cross_product
    DOUBLE PRECISION                            :: cross1, cross2, cross3

    END FUNCTION cross_product
    end interface

    !Input Variables
    type(Vortex_Segment), DIMENSION (:), INTENT(INOUT)  ::  VtxSgmntSet ! Corresponding Vtx Segments set
    type(Nodal_Point), DIMENSION (:), INTENT(IN)     ::  NodesSet    ! Corresponding Nodes set
    type(Panel), DIMENSION (:), INTENT(INOUT)        ::  PnlSet      ! Corresponding Panels set

    !Local Variables
    DOUBLE PRECISION, DIMENSION(3) :: N1    ! Node coords extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: N2    ! Node coords extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: N3    ! Node coords extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: N4    ! Node coords extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: VS1   ! Segment extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: VS2   ! Segment extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: VS3   ! Segment extraction variable
    DOUBLE PRECISION, DIMENSION(3) :: VS4   ! Segment extraction variable
    INTEGER            :: nPanels  ! Number of panels of the set
    DOUBLE PRECISION               :: ShArea !Shared Area aux variable
    INTEGER            :: i,j,k    ! Loop index

    nPanels = size(PnlSet)

    DO i=1,nPanels

        !Panel AREA
        select case (PnlSet(i)%Panel_Type)

            case default
            !Nodes Extraction
            N1 = NodesSet(PnlSet(i)%PanelNodes(1))%xyz
            N2 = NodesSet(PnlSet(i)%PanelNodes(2))%xyz
            N3 = NodesSet(PnlSet(i)%PanelNodes(3))%xyz
            N4 = NodesSet(PnlSet(i)%PanelNodes(4))%xyz

            !Area Segments calculation
            VS1 = N2 - N1
            VS2 = N3 - N2
            VS3 = N3 - N4
            VS4 = N4 - N1

            PnlSet(i)%Area = 0.5*(norm2(cross_product(VS1,VS2)) + norm2(cross_product(VS3,VS4)))
            !print*,'Area = ', Area

        case (3)
            !Nodes Extraction
            N1 = NodesSet(PnlSet(i)%PanelNodes(1))%xyz
            N2 = NodesSet(PnlSet(i)%PanelNodes(2))%xyz
            N3 = NodesSet(PnlSet(i)%PanelNodes(3))%xyz

            !Area Segments calculation
            VS1 = N2 - N1
            VS2 = N3 - N1

            PnlSet(i)%Area = 0.5*(norm2(cross_product(VS1,VS2)))
            !print*,'Area = ', Area
        end select

        !PANEL VtxSEGMENTS SHARED AREA
        DO k=1,PnlSet(i)%nPanelSegments
            ShArea = 0.0
            DO j=1,VtxSgmntSet(abs(PnlSet(i)%PanelSegments(k)))%nSharingPanels
                ShArea= ShArea + PnlSet(abs(VtxSgmntSet(abs(PnlSet(i)%PanelSegments(k)))%SharingPanels(j)))%Area
                !print *, 'ShArea:  ', ShArea
            END DO
            VtxSgmntSet(abs(PnlSet(i)%PanelSegments(k)))%SharedArea = ShArea
        END DO
    END DO


    END SUBROUTINE Area_calc

    END MODULE classPanel