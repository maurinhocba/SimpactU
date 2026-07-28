    MODULE AeroMatrix

    ! Gathers subroutines for [A] Matrix and [RHS] vector for grid/wake input, includes computation and printing procedures.


    USE classVortex_Segment
    USE classNodal_Point
    USE classPanel
    USE P_S_N
    USE classGrid
    USE classWake
    USE classBody


    IMPLICIT NONE
    PUBLIC

    CONTAINS ! ==================================================================

    SUBROUTINE A_Mat_Assembler (BodySet, nBodies, Full_AMat)
    ! Calculates the Influence Coeficient Matrix for a pair of Grids (using 'InfCoefMatrix', APPLIES FOR GRIDS AND BC_GRIDS)
    ! and assambles it into the Full_A_Matrix, according to the AMatLoc index of each grid.

    !Input Variables
    type(BODY), ALLOCATABLE, DIMENSION(:), INTENT(IN)  :: BodySet ! Gives the set of grids for the matrix constructor.
    INTEGER                                            :: nBodies  ! Number of bodies of the simulation.

    !Output Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: Full_AMat    ! Global Influence Coeficient Matrix.

    !Local Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A_Matrix
    INTEGER                           :: i,j,k,l ! loop index
    INTEGER                           :: iLoc, jLoc, iStride, jStride ! Position index extraction variables

    
    !Init
    Full_AMat = 0.d0
    
    ! Grids
    !print *,'Bodies first loop ',BodySet(1)%nGrids, BodySet(1)%nBCGrids
    DO i=1,nBodies !  First loop on Bodies
        DO j=1,BodySet(i)%nGrids !  First loop on Grids
            iLoc = BodySet(i)%Grids(j)%AMatLoc ! i position extraction
            iStride = BodySet(i)%Grids(j)%nPanels - 1 ! i Stride calculation
            DO k=1,nBodies !  Second loop on Bodies
                DO l=1,BodySet(k)%nGrids !  Second loop on Grids
                    !iLoc = BodySet(i)%Grids(j)%AMatLoc ! i position extraction
                    !iStride = BodySet(i)%Grids(j)%nPanels - 1 ! i Stride calculation
                    jLoc = BodySet(k)%Grids(l)%AMatLoc ! j position extraction
                    jStride = BodySet(k)%Grids(l)%nPanels - 1 ! j Stride calculation
                    call InfCoefMatrix (BodySet(i)%Grids(j),BodySet(k)%Grids(l),A_Matrix) ! Grid on Grid A_Matrix calculation (Local Matrix)
                    !Full_AMat(iLoc:iLoc+iStride,jLoc:jLoc+jStride) = A_Matrix  ! Local matrix assembly
                    Full_AMat(jLoc:jLoc+jStride,iLoc:iLoc+iStride) = A_Matrix  ! Local matrix assembly
                    deallocate(A_Matrix)  ! Local Matrix deallocation
                END DO
                !print *, 'BCGrids loop'
                DO l=1,BodySet(k)%nBCGrids !  Second loop on BCGrids
                    !iLoc = BodySet(i)%BCGrids(j)%AMatLoc ! i position extraction
                    !iStride = BodySet(i)%BCGrids(j)%nPanels - 1 ! i Stride calculation
                    jLoc = BodySet(k)%BCGrids(l)%AMatLoc ! j position extraction
                    jStride = BodySet(k)%BCGrids(l)%nPanels - 1 ! j Stride calculation
                    call InfCoefMatrix (BodySet(i)%Grids(j),BodySet(k)%BCGrids(l),A_Matrix) ! Grid on Grid A_Matrix calculation (Local Matrix)
                    !Full_AMat(iLoc:iLoc+iStride,jLoc:jLoc+jStride) = A_Matrix  ! Local matrix assembly
                    Full_AMat(jLoc:jLoc+jStride,iLoc:iLoc+iStride) = A_Matrix  ! Local matrix assembly
                    deallocate(A_Matrix)  ! Local Matrix deallocation
                    !print *,'deallocated'
                END DO
            END DO
        END DO
        DO j=1,BodySet(i)%nBCGrids !  First loop on BCGrids
            iLoc = BodySet(i)%BCGrids(j)%AMatLoc ! i position extraction
            iStride = BodySet(i)%BCGrids(j)%nPanels - 1 ! i Stride calculation
            DO k=1,nBodies !  Second loop on Bodies
                DO l=1,BodySet(k)%nGrids !  Second loop on Grids and BCGrids
                    jLoc = BodySet(k)%Grids(l)%AMatLoc ! j position extraction
                    jStride = BodySet(k)%Grids(l)%nPanels - 1 ! j Stride calculation
                    call InfCoefMatrix (BodySet(i)%BCGrids(j),BodySet(k)%Grids(l),A_Matrix) ! Grid on Grid A_Matrix calculation (Local Matrix)
                    !Full_AMat(iLoc:iLoc+iStride,jLoc:jLoc+jStride) = A_Matrix  ! Local matrix assembly
                    Full_AMat(jLoc:jLoc+jStride,iLoc:iLoc+iStride) = A_Matrix  ! Local matrix assembly
                    deallocate(A_Matrix)  ! Local Matrix deallocation
                    !print *,'deallocated'
                END DO
                DO l=1,BodySet(k)%nBCGrids !  Second loop on Grids
                    !iLoc = BodySet(i)%BCGrids(j)%AMatLoc ! i position extraction
                    jLoc = BodySet(k)%BCGrids(l)%AMatLoc ! j position extraction
                    !iStride = BodySet(i)%BCGrids(j)%nPanels - 1 ! i Stride calculation
                    jStride = BodySet(k)%BCGrids(l)%nPanels - 1 ! j Stride calculation
                    call InfCoefMatrix (BodySet(i)%BCGrids(j),BodySet(k)%BCGrids(l),A_Matrix) ! Grid on Grid A_Matrix calculation (Local Matrix)
                    !Full_AMat(iLoc:iLoc+iStride,jLoc:jLoc+jStride) = A_Matrix  ! Local matrix assembly
                    Full_AMat(jLoc:jLoc+jStride,iLoc:iLoc+iStride) = A_Matrix  ! Local matrix assembly
                    deallocate(A_Matrix)  ! Local Matrix deallocation
                END DO
            END DO
        END DO
    END DO

    !! BC_Grids
    !print *, 'BC_grids AMat:---'
    !DO i=1,nBodies !  First loop on Bodies
    !    print *, 'i, nBCgrids  : ',i,BodySet(i)%nBCGrids
    !    DO j=1,BodySet(i)%nBCGrids !  First loop on Grids
    !        DO k=1,nBodies !  Second loop on Bodies
    !            DO l=1,BodySet(k)%nBCGrids !  Second loop on Grids
    !                iLoc = BodySet(i)%BCGrids(j)%AMatLoc ! i position extraction
    !                jLoc = BodySet(k)%BCGrids(l)%AMatLoc ! j position extraction
    !                iStride = BodySet(i)%BCGrids(j)%nPanels - 1 ! i Stride calculation
    !                jStride = BodySet(k)%BCGrids(l)%nPanels - 1 ! j Stride calculation
    !                call InfCoefMatrix (BodySet(i)%BCGrids(j),BodySet(k)%BCGrids(l),A_Matrix) ! Grid on Grid A_Matrix calculation (Local Matrix)
    !                Full_AMat(iLoc:iLoc+iStride,jLoc:jLoc+jStride) = A_Matrix  ! Local matrix assembly
    !                deallocate(A_Matrix)  ! Local Matrix deallocation
    !            END DO
    !        END DO
    !    END DO
    !END DO

    !print *, 'FULL_A_Mat---: ',Full_AMat(:,:)
    !pause

    END SUBROUTINE A_Mat_Assembler

    SUBROUTINE AeroMat_GenLoc (BodySet, nBodies, Full_AMat, Full_RHS)
    ! Sorts Grids and BC_Grids of the simulation so as to give each one a position into the global A_Matrix.
    ! This position is represented by the AMatLoc variable of the grid which gives the index of the first
    ! element for each matrix inside the global matrix.
    ! The position given is held until the end of the simulation and applies also for the RHS vector.
    ! For convenience Grids are positioned before BC_Grids, following bodies numeration.
    ! The same is applied for the Full_RHS vector, for positions and allocation
    ! Full_AMat and Full_RHS are generated and allocated.

    !Input Variables
    INTEGER, INTENT(IN) :: nBodies !Number of Bodies of the simulation

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(INOUT) :: BodySet ! Set of bodies of the simulation.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: Full_AMat    ! Global Influence Coeficient Matrix.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: Full_RHS       ! Global Right Hand Side Vector.

    !Local Variables
    INTEGER :: i,j ! Loop Index
    INTEGER :: Posic ! Aux positioning variable

    Posic = 1 ! Set counter to 1, first index of matrix. Actualization is cummulative over bodies.
    print *,'Grids loop.....'
    DO i=1,nBodies  ! Bodies loop
        DO j=1,BodySet(i)%nGrids !Grids loop
            BodySet(i)%Grids(j)%AMatLoc = Posic ! Grid positioning
            Posic = Posic + BodySet(i)%Grids(j)%nPanels ! The index stride is nPanels of the previous grid.
        END DO
        print *,'BCGrids loop.....'
        DO j=1,BodySet(i)%nBCGrids ! BC_Grids loop
            BodySet(i)%BCGrids(j)%AMatLoc = Posic ! Grid positioning
            Posic = Posic + BodySet(i)%BCGrids(j)%nPanels ! The index stride is nPanels of the previous grid.
        END DO
    END DO
    print *,'ALLOCATION A MAT', Posic
    ALLOCATE (Full_AMat (Posic-1,Posic-1))
    ALLOCATE (Full_RHS (Posic-1))
    print *,'end AeroMat_GenLoc'
    END SUBROUTINE AeroMat_GenLoc


    SUBROUTINE InfCoefMatrix (Grd1, Grd2, A_Matrix)
    ! USES POLYMORPHIC VARIABLES FOR Grd1, Grd2 SO AS TO RECIEVE ALSO BC_Grids AND EXECUTE THE SAME PROCESSING.
    !Calculates the Influence Coeficient Matrix (A_Matrix) of Grid1(emitting) over Grid2(receiving),
    !defined as the array with the normal-wise speed induced in each control point by each panel vortex ring, with Gamma=1.
    !If Grid1=Grid2 the A_Matrix represents the influence of the grid on itself.
    !Gives a block to assemble into the Full A_Matrix with the correspondig subroutine.


    interface
    FUNCTION VtxSgmnt_IndSpeed (NP1, NP2, CP, Gamma)
    !Calculates the speed induced by a vortex segment (from NP1 to NP2) of intensity Gamma, in the control point CP.

    IMPLICIT NONE

    !Input Variables
    DOUBLE PRECISION, DIMENSION (3), INTENT(IN)  :: NP1, NP2    ! Nodal Points of the Segment, first-last.
    DOUBLE PRECISION, DIMENSION (3), INTENT(IN)  :: CP          ! Control Point to calculate induced speed in.
    DOUBLE PRECISION, INTENT(IN)                 :: Gamma       ! Circulation of the segment.

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
    CLASS(GRID), INTENT(IN) :: Grd1 ! Emitting, induces the speed over Grd2 CtrlPoints. APPLIES FOR GRIDS AND BC_GRIDS
    CLASS(GRID), INTENT(IN) :: Grd2 ! Receiving, gives the CtrlPoints for induced speed calculation from Grd1 vorticity. APPLIES FOR GRIDS AND BC_GRIDS
    !type(GRID), INTENT(IN) :: Grd1 ! Emitting, induces the speed over Grd2 CtrlPoints.
    !type(GRID), INTENT(IN) :: Grd2 ! Receiving, gives the CtrlPoints for induced speed calculation from Grd1 vorticity.

    !Output Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: A_Matrix

    !Local Variables
    INTEGER                                        :: NumberOfPanels1,NumberOfPanels2  ! Number of panels of the grids. Variable extraction.
    !INTEGER                            :: PnlType         ! Panel Type code. Variable association.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Xa           ! Panel Nodes Coordinates Matrix
    INTEGER                                        :: NodeId       ! Node Label. Variable extraction.
    DOUBLE PRECISION, DIMENSION (3)                :: CPCoor       ! Coordinates of Control Point.
    DOUBLE PRECISION, DIMENSION (3)                :: CPNormal     ! Normal Versor of the panel at Control Point.
    DOUBLE PRECISION, DIMENSION (3)                :: IndSpeed     ! Speed induced at control point by each panel circulation.
    DOUBLE PRECISION, DIMENSION (3)                :: CPIndSpeed   ! Cumulative induced speed at Control Point, over every panel.
    INTEGER :: i,j,n           ! Loops indexes



    !Matrix Allocation
    NumberOfPanels1 = Grd1%nPanels
    !print *,'G1 nPanels = ',NumberOfPanels1
    NumberOfPanels2 = Grd2%nPanels
    !print *,'G2 nPanels = ',NumberOfPanels2
    Allocate(A_Matrix(NumberOfPanels2,NumberOfPanels1))
               
    DO j=1,NumberOfPanels1      !Panel loop. Emitting
        !PType: associate (PnlType => Grd1%PanelSet(j)%Panel_Type) ! Panel Type code. Variable association
            !nPnlNodes: associate (nNodes => Grd1%PanelSet(j)%nPanelNodes)

                !CPIndSpeed = 0.d0 !Sets sum to zero.

                call InGrid_Xa (Grd1%PanelSet(j)%PanelNodes, Grd1%NodalPointSet, Grd1%PanelSet(j)%Panel_Type, Xa) ! Xa function call.
                
				 !print *,'Grd2%PanelSet(j)%Panel_Type = ', Grd2%PanelSet(j)%Panel_Type 
				
                DO i=1,NumberOfPanels2    ! CP loop. Receiving
                    CPCoor = Grd2%PanelSet(i)%CtrlPoint%xyz
                    CPNormal = Grd2%PanelSet(i)%CtrlPoint%NormalVersor
                    !select case (PnlType)

                    !case (1) ! 4-node Cuad

                        IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,4), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,1), CPCoor, 1.d0)
                        !Sum of the induced speed by each vtx sgmnt at the control point.

                    !case (3) ! 3-node Triag
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,1), CPCoor, 1.d0)
                    !
                    !case (8) ! 8-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,5), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,5), Xa(:,2), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,6), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,6), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,7), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,7), Xa(:,4), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,8), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,8), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (51) ! 5-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,5), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,5), Xa(:,2), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,4), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (52) ! 5-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,6), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,6), Xa(:,3), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,4), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !case (53) ! 5-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,7), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,7), Xa(:,4), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !case (54) ! 5-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,4), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,8), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,8), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (61) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,5), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,5), Xa(:,2), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,6), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,6), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,4), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (62) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,6), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,6), Xa(:,3), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,7), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,7), Xa(:,4), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (63) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,7), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,7), Xa(:,4), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,8), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,8), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (64) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,5), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,5), Xa(:,2), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,4), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,8), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,8), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (71) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,5), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,5), Xa(:,2), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,6), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,6), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,7), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,7), Xa(:,4), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (72) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,2), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,6), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,6), Xa(:,3), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,7), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,7), Xa(:,4), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,8), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,8), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (73) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,5), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,5), Xa(:,2), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,3), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,7), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,7), Xa(:,4), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,8), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,8), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.
                    !
                    !case (74) ! 6-node Cuad
                    !
                    !    IndSpeed=VtxSgmnt_IndSpeed (Xa(:,1), Xa(:,5), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,5), Xa(:,2), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,2), Xa(:,6), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,6), Xa(:,3), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,3), Xa(:,4), CPCoor, 1.d0) + VtxSgmnt_IndSpeed (Xa(:,4), Xa(:,8), CPCoor, 1.d0) + &
                    !        VtxSgmnt_IndSpeed (Xa(:,8), Xa(:,1), CPCoor, 1.d0)
                    !    !Sum of the induced speed by each vtx sgmnt at the control point.

                    !end select
                      
                    A_Matrix(i,j)=dot_product(IndSpeed,CPNormal)

                    !CPIndSpeed = CPIndSpeed + IndSpeed
                END DO
                     
            !end associate nPnlNodes

        !end associate PType

        deallocate (Xa)
        !Grd%PanelSet(i)%CtrlPoint%InducedVelocity = CPIndSpeed   !CHECK!! must be actualized at every step

    END DO


    END SUBROUTINE InfCoefMatrix


    SUBROUTINE InfCoefMat_Printing (A_Matrix, FileUnit)
    !Prints the Influence Coeficien Matrix [A] in a file specified by FileUnit.

    !Input Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: A_Matrix
    INTEGER, INTENT (IN)                          :: FileUnit

    !Local Variables
    INTEGER :: rows, clmns ! Rows and Columns of the matrix.
    INTEGER :: i,j         ! Loop indexes.


    rows = size(A_Matrix,1)
    clmns = size(A_Matrix,2)

    write (FileUnit,'(A11,/)') 'A_Matrix: '
    write (FileUnit,'(A21,I7,A1,I7)') 'Size (rowsXcolumns): ',rows,'x',clmns

    DO i=1,rows
        DO j=1,clmns
            write (FileUnit, '(F20.16)',advance='no') A_Matrix(i,j)
        END DO
        write (FileUnit, '(/)')
    END DO

    END SUBROUTINE InfCoefMat_Printing
    
    SUBROUTINE RHSVector(Grd, BodySet, V_inf, RHS)
    !Calculates the Right Hand Side vector of a Grid, defined as the array with the normal-wise speed at each
    !control point produced by freestream, body rotation/traslation speed and wakes influence. Influence of the grid
    !on itself is not considered.
	!Also stores the wakes influence induced speed at each Control Point.
    ! USES POLYMORPHIC VARIABLES SO AS TO PROCESS GRIDS AND BC_GRIDS.

    ! The resulting RHS must be assembled into the global or Full vector.

    !REQUIRES THE COMPUTATION OF WIND VELOCITY AND OWN VELOCITY AT EACH CONTROL POINT -> CONTROL POINT SUBROUTINE
    !CP_Velocity -> Subroutine for storing velocities at control point.

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
    DOUBLE PRECISION, PARAMETER      :: PI=3.1415926535897933d0
    DOUBLE PRECISION                 :: L               ! Length of the vortex segment.
    DOUBLE PRECISION                 :: delta    ! "Cutoff" radius parameter (VAN GARREL)
    END FUNCTION VtxSgmnt_IndSpeed
    end interface

    !Input Variables
    !class(GRID), INTENT(IN)                           :: Grd !Grid over which the RHS is calculated. APPLIES FOR GRIDs and BC_GRIDS
    type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(IN) :: BodySet ! Set of bodies of the simulation. Gives all Wakes_Sets to use
    !type(WAKE), INTENT (IN) :: Wke
    DOUBLE PRECISION,  DIMENSION (3), INTENT(IN)     :: V_Inf
    
    !Output Variables
    CLASS(GRID), INTENT(INOUT)                               :: Grd ! Grid over which the RHS is calculated. The CtrlPoints induced speed by wakes is also stored. APPLIES FOR GRIDs and BC_GRIDS
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: RHS

    !Local Variables
    INTEGER             :: NumberOfPanels   !Number of Panels of the grid.
    INTEGER             :: nBodies ! Number of Bodies of the simulation.
    INTEGER             :: i,j,k  !loop indexes
    DOUBLE PRECISION, DIMENSION (3) :: IndSpd ! Speed induced by one Wake Segment at a CP.
    DOUBLE PRECISION, DIMENSION (3) :: WkeIndSpeed, SumWkeIndSpeed ! Induced Speed by wake at a CP.
    
    DOUBLE PRECISION, DIMENSION (3) :: CPxyz, OwnV, WindV, Norm ! Aux Variables, Replaces Associate construct
    
    !INTEGER             :: N1, N2 ! Segment nodes ID

    nBodies = size(BodySet) ! Number of Bodies extraction.

    !Vector Allocation
    NumberOfPanels = Grd%nPanels
    Allocate(RHS(NumberOfPanels))
      
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(RHS, GRD, BodySet, nBodies, NumberOfPanels, V_inf)
    !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(CPxyz, OwnV, WindV, Norm)

    DO i=1,NumberOfPanels
        SumWkeIndSpeed=[0.d0,0.d0,0.d0]

        CPxyz = Grd%PanelSet(i)%CtrlPoint%xyz
        OwnV  = Grd%PanelSet(i)%CtrlPoint%OwnVelocity 
        !WindV = Grd%PanelSet(i)%CtrlPoint%WindVelocity
        WindV = V_inf
        Norm  = Grd%PanelSet(i)%CtrlPoint%NormalVersor
        
        !print *,'OwnV@CP >>: ',OwnV
                        !pause
                        DO j=1,nBodies
                            DO k=1,BodySet(j)%nWakes
                                call Wake_Influence_IndSpeed (BodySet(j)%Wakes(k), CPxyz, WkeIndSpeed)
                                SumWkeIndSpeed = SumWkeIndSpeed + WkeIndSpeed
                                !print *, 'CPxyz, k: ', CPxyz, k
                                !print *, 'SUMWke indspeed: ', SumWkeIndSpeed
                            END DO
						END DO

						Grd%PanelSet(i)%CtrlPoint%WKEInducedVelocity = SumWkeIndSpeed  ! Stores the	velocity that wakes induce at the control point to avoid recalculation.
						
                        select type (Grd) ! Type selector: consider the case of processing a BC_GRID in which the transporation velocity must be considered in the RHS.
                        type is (BC_Grid)
                            !print *, 'BC type windV: ', WindV
                            !pause
                            if (Grd%TrnspCondition==1) then  ! Case for Transp velocity constant (non-zero) over the grid.
                                RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm) + Grd%TranspirationVel
                            else if (Grd%TrnspCondition==2) then ! Case for Transp velocity variable for each panel.
                                RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm) + Grd%PanelSet(i)%CtrlPoint%TrnspVelocity
                            else ! Case BC_GRID with no penetration condition
                                RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm)
                                !print *, 'Wind: ', WindV
                                !print *, 'OwnV: ', OwnV
                                !print *, 'Norm: ', Norm
                                !print *, 'SUMWke indspeed: ', SumWkeIndSpeed
                                !print *, 'RHSi: ', RHS(i)
                            end if
                            class default  ! Default case: type (GRID)
                            !print *, 'Grid type SumWkeIndSpeed: ', SumWkeIndSpeed
                            !pause
                            RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm)
                        end select


    END DO

    !$OMP END DO
    !$OMP END PARALLEL

	END SUBROUTINE RHSVector

    !SUBROUTINE RHSVector_old (Grd, BodySet, RHS)
    !!Calculates the Right Hand Side vector of a Grid, defined as the array with the normal-wise speed at each
    !!control point produced by freestream, body rotation/traslation speed and wakes influence. Influence of the grid
    !!on itself is not considered.
    !! USES POLYMORPHIC VARIABLES SO AS TO PROCESS GRIDS AND BC_GRIDS.
    !
    !! The resulting RHS must be assembled into the global or Full vector.
    !
    !!REQUIRES THE COMPUTATION OF WIND VELOCITY AND OWN VELOCITY AT EACH CONTROL POINT -> CONTROL POINT SUBROUTINE
    !!CP_Velocity -> Subroutine for storing velocities at control point.
    !
    !interface
    !FUNCTION VtxSgmnt_IndSpeed (NP1, NP2, CP, Gamma)
    !!Calculates the speed induced by a vortex segment (from NP1 to NP2) of intensity Gamma, in the control point CP.
    !
    !IMPLICIT NONE
    !
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
    !!Input Variables
    !CLASS(GRID), INTENT(IN)                          :: Grd ! Grid over which the RHS is calculated. APPLIES FOR GRIDs and BC_GRIDS
    !!type(GRID), INTENT(IN)                           :: Grd ! Grid over which the RHS is calculated
    !type(BODY), ALLOCATABLE, DIMENSION(:),INTENT(IN) :: BodySet ! Set of bodies of the simulation. Gives all Wakes_Sets to use
    !!type(WAKE), INTENT (IN) :: Wke
    !
    !!Output Variables
    !DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: RHS
    !
    !!Local Variables
    !INTEGER             :: NumberOfPanels   !Number of Panels of the grid.
    !INTEGER             :: nBodies ! Number of Bodies of the simulation.
    !INTEGER             :: i,j,k  !loop indexes
    !DOUBLE PRECISION, DIMENSION (3) :: IndSpd ! Speed induced by one Wake Segment at a CP.
    !DOUBLE PRECISION, DIMENSION (3) :: WkeIndSpeed, SumWkeIndSpeed ! Induced Speed by wake at a CP.
    !!INTEGER             :: N1, N2 ! Segment nodes ID
    !
    !nBodies = size(BodySet) ! Number of Bodies extraction.
    !
    !!Vector Allocation
    !NumberOfPanels = Grd%nPanels
    !Allocate(RHS(NumberOfPanels))
    !
    !DO i=1,NumberOfPanels
    !    SumWkeIndSpeed=[0.0,0.0,0.0]
    !
    !    !Own and Wind Speed influence
    !    Coords: associate (CPxyz => Grd%PanelSet(i)%CtrlPoint%xyz) ! Variable association
    !        Own: associate (OwnV => Grd%PanelSet(i)%CtrlPoint%OwnVelocity) ! Variable association
    !            Wind: associate (WindV => Grd%PanelSet(i)%CtrlPoint%WindVelocity) ! Variable association
    !                Normal: associate (Norm => Grd%PanelSet(i)%CtrlPoint%NormalVersor) ! Variable association
    !
    !                    !print *,'OwnV@CP >>: ',OwnV
    !                    !pause
    !                    DO j=1,nBodies
    !                        DO k=1,BodySet(j)%nWakes
    !                            call Wake_Influence_IndSpeed (BodySet(j)%Wakes(k), CPxyz, WkeIndSpeed)
    !                            SumWkeIndSpeed = SumWkeIndSpeed + WkeIndSpeed
    !                            !print *, 'CPxyz, k: ', CPxyz, k
    !                            !print *, 'SUMWke indspeed: ', SumWkeIndSpeed
    !                        END DO
    !                    END DO
    !
    !                    select type (Grd) ! Type selector: consider the case of processing a BC_GRID in which the transporation velocity must be considered in the RHS.
    !                    type is (BC_Grid)
    !                        !print *, 'BC type windV: ', WindV
    !                        !pause
    !                        if (Grd%TrnspCondition==1) then  ! Case for Transp velocity constant (non-zero) over the grid.
    !                            RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm) + Grd%TranspirationVel
    !                        else if (Grd%TrnspCondition==2) then ! Case for Transp velocity variable for each panel.
    !                            RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm) + Grd%PanelSet(i)%CtrlPoint%TrnspVelocity
    !                        else ! Case BC_GRID with no penetration condition
    !                            RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm)
    !                            !print *, 'Wind: ', WindV
    !                            !print *, 'OwnV: ', OwnV
    !                            !print *, 'Norm: ', Norm
    !                            !print *, 'SUMWke indspeed: ', SumWkeIndSpeed
    !                            !print *, 'RHSi: ', RHS(i)
    !                        end if
    !                        class default  ! Default case: type (GRID)
    !                        !print *, 'Grid type SumWkeIndSpeed: ', SumWkeIndSpeed
    !                        !pause
    !                        RHS(i)= -dot_product((WindV - OwnV + SumWkeIndSpeed),Norm)
    !                    end select
    !
    !                end associate Normal
    !            end associate Wind
    !        end associate Own
    !    end associate Coords
    !
    !
    !END DO
    !
    !END SUBROUTINE RHSVector_old

    SUBROUTINE RHS_Assembler (BodySet, nBodies, V_inf, Full_RHS)
    ! Calculates the RHS_Vector for a Grid (using 'RHSVector', APPLIES FOR GRIDS AND BC_GRIDS)
    ! and assambles it into the Full_RHS, according to the AMatLoc index of each grid.

    !Input Variables
    INTEGER                                            :: nBodies  ! Number of bodies of the simulation.
    DOUBLE PRECISION,  DIMENSION (3), INTENT(IN)       :: V_Inf

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: BodySet ! Gives the set of grids for the matrix constructor.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: Full_RHS    ! Global RHS Vector.

    !Local Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: RHS ! One Grid RHS Vector
    INTEGER                           :: i,j,k,l ! loop index
    INTEGER                           :: iLoc, iStride ! Position index extraction variables

    ! Grids
    DO i=1,nBodies !  Loop on Bodies
        DO j=1,BodySet(i)%nGrids !  Loop on Grids
            iLoc = BodySet(i)%Grids(j)%AMatLoc ! i position extraction
            iStride = BodySet(i)%Grids(j)%nPanels - 1 ! i Stride calculation
            call RHSVector (BodySet(i)%Grids(j),BodySet,V_inf,RHS) ! Single Grid RHS vector calculation ("Local" vector)
            !print *,'RHS vect : ', RHS(:)
            Full_RHS(iLoc:iLoc+iStride) = RHS  ! Local Vector assembly
            deallocate(RHS)  ! Local Vector deallocation
        END DO
    END DO

    ! BC_Grids
    DO i=1,nBodies !  Loop on Bodies
        DO j=1,BodySet(i)%nBCGrids !  Loop on Grids
            iLoc = BodySet(i)%BCGrids(j)%AMatLoc ! i position extraction
            iStride = BodySet(i)%BCGrids(j)%nPanels - 1 ! i Stride calculation
            call RHSVector (BodySet(i)%BCGrids(j),BodySet,V_inf,RHS) ! Single Grid RHS vector calculation ("Local" vector)
            Full_RHS(iLoc:iLoc+iStride) = RHS  ! Local Vector assembly
            deallocate(RHS)  ! Local Vector deallocation
        END DO
    END DO


    END SUBROUTINE RHS_Assembler

    SUBROUTINE RHSv_Printing (RHS, FileUnit)
    !Prints the Right Hand Side Vector [RHS] in a file specified by FileUnit.

    !Input Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: RHS
    INTEGER, INTENT (IN)                        :: FileUnit

    !Local Variables
    INTEGER  :: rows  ! Rows of the vector.
    INTEGER  :: i     ! Loop index.

    rows = size(RHS)

    write (FileUnit,'(A11,/)') 'RHS vector: '
    write (FileUnit,'(A6,I7)') 'Size: ',rows

    DO i=1,rows
        write (FileUnit, '(F20.16)') RHS(i)
    END DO



    END SUBROUTINE RHSv_Printing

    SUBROUTINE Vorticity_Process (BodySet, nBodies, G_vect)
    !Transfers Vorticity obtained from the Solver (G_vect) and sets it into Panels and Segments of Grids and BC_Grids
    !Uses AMatLoc index for identifing G_vect positions in accordance with grids.
    !Calls 'PanelSet_RingVtct_Loading' and 'RingVtct_2_SgmntCirc' (Grids' subroutines)

    !Input Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:),INTENT(IN)   :: G_vect  !Global vorticity vector obtained from Solver.
    INTEGER, INTENT(IN)                          :: nBodies !Number of Bodies of the simulation

    !Output Variables
    type(BODY), ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: BodySet ! Gives the set of grids and bc_grids.

    !Local Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: G ! Local G_vector for each grid
    INTEGER      :: i,j,k,l ! loop index
    INTEGER      :: iLoc, iStride ! Position index extraction variables

    ! Grids
    DO i=1,nBodies !  Loop on Bodies
        DO j=1,BodySet(i)%nGrids !  Loop on Grids
            iLoc = BodySet(i)%Grids(j)%AMatLoc ! i position extraction
            iStride = BodySet(i)%Grids(j)%nPanels - 1 ! i Stride calculation
            allocate (G (BodySet(i)%Grids(j)%nPanels)) ! Local Vector allocation
            G = G_vect(iLoc:iLoc+iStride)   ! Local G vector extraction
            call PanelSet_RingVtct_Loading (BodySet(i)%Grids(j), G)
            call RingVrtct_2_SgmntCirc (BodySet(i)%Grids(j))
            deallocate(G)  ! Local Vector deallocation
        END DO
    END DO

    ! BCGrids
    DO i=1,nBodies !  Loop on Bodies
        DO j=1,BodySet(i)%nBCGrids !  Loop on Grids
            iLoc = BodySet(i)%BCGrids(j)%AMatLoc ! i position extraction
            iStride = BodySet(i)%BCGrids(j)%nPanels - 1 ! i Stride calculation
            allocate (G (BodySet(i)%BCGrids(j)%nPanels)) ! Local Vector allocation
            G = G_vect(iLoc:iLoc+iStride)   ! Local G vector extraction
            call PanelSet_RingVtct_Loading (BodySet(i)%BCGrids(j), G)
            call RingVrtct_2_SgmntCirc (BodySet(i)%BCGrids(j))
            deallocate(G)  ! Local Vector deallocation
        END DO
    END DO

    END SUBROUTINE Vorticity_Process

    SUBROUTINE Gv_Printing (G, FileUnit)
    !Prints the Right Hand Side Vector [RHS] in a file specified by FileUnit.

    !Input Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(IN) :: G
    INTEGER, INTENT (IN)                        :: FileUnit

    !Local Variables
    INTEGER  :: rows  ! Rows of the vector.
    INTEGER  :: i     ! Loop index.

    rows = size(G)

    write (FileUnit,'(A11,/)') 'G vector: '
    write (FileUnit,'(A6,I7)') 'Size: ',rows

    DO i=1,rows
        write (FileUnit, '(F25.16)') G(i)
    END DO

    END SUBROUTINE Gv_Printing

    END MODULE