    MODULE classControl_Point

    ! "Control_Point" extended class type definition
    ! Specialization of "Point"

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    USE classPoint, only: Point
    IMPLICIT NONE
    PUBLIC

    interface
    FUNCTION cross_product(U,V)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
    DOUBLE PRECISION, DIMENSION(3)              :: cross_product
    DOUBLE PRECISION                            :: cross1, cross2, cross3

    END FUNCTION cross_product
    end interface



    TYPE, EXTENDS(Point), Public ::    Control_Point ! extended type of "Point"
        
    !Attributes:
    DOUBLE PRECISION, DIMENSION (3)     :: OwnVelocity         ! control point own velocity
    DOUBLE PRECISION, DIMENSION (3)     :: InducedVelocity     ! TOTAL induced velocity at control point 
	DOUBLE PRECISION, DIMENSION (3)     :: WKEInducedVelocity  ! induced velocity at control point BY WAKES ONLY
    DOUBLE PRECISION, DIMENSION (3)     :: WindVelocity        ! free stream velocity at control point
    DOUBLE PRECISION, DIMENSION (3)     :: NormalVersor        ! Panel outgoing normal versor at cotrol point
    DOUBLE PRECISION                    :: TrnspVelocity       ! Transpiration velocity. Always normal to surface, is considered the proyection over Normal Versor.
    DOUBLE PRECISION                    :: Delta_Cp            ! Pressure Coefficient Difference across the lifting surface at the control point.

    CONTAINS ! methods' names

    PROCEDURE, Public, NoPass :: CP_Coords
    PROCEDURE, Public, NoPass :: CP_NormalVersor
    PROCEDURE, Public, NoPass :: CP_Velocity

    END TYPE Control_Point

    CONTAINS ! ==================================================================

    SUBROUTINE CP_Coords (Xa, PShapeFunct, CP_xyz)
    ! Calculates the CtrlPoint coordinates for a panel

    !Input Variables
    DOUBLE PRECISION, DIMENSION (:,:), INTENT(IN) :: Xa             ! Panel Nodes Coords
    DOUBLE PRECISION, DIMENSION (8), INTENT(IN)   :: PShapeFunct    ! Panel Shape Functions

    !Output Variables
    DOUBLE PRECISION, DIMENSION (3), INTENT(OUT)  :: CP_xyz         ! Output vector of coords

    !Local Variables
    INTEGER ::  i ! Loop index

    CP_xyz = 0 ! Initialize with zeros

    do i=1,8

        if (PShapeFunct(i) /= 0) then               ! Shape Functions vector is fill with zeros according to panel type
            CP_xyz = CP_xyz + Xa(:,i)*PShapeFunct(i)    ! Undefined size vector dot product: X_cp=[Xa]*[N]_0
        end if

    end do

    END SUBROUTINE CP_Coords

    SUBROUTINE CP_NormalVersor (Xa, Derived_PSF, NormalVersor)

    !Input Variables
    DOUBLE PRECISION, DIMENSION (:,:), INTENT(IN)   :: Xa             ! Panel Nodes Coords
    DOUBLE PRECISION, DIMENSION (8,2), INTENT(IN)   :: Derived_PSF    ! Panel Shape Functions

    !Output Variables
    DOUBLE PRECISION, DIMENSION (3), INTENT(OUT)    :: NormalVersor   ! Output Normal Versor

    !Local Variables
    DOUBLE PRECISION, DIMENSION (3,2)               :: JacM           ! Jacobian Matrix of the Panel Transform
    DOUBLE PRECISION, DIMENSION (3)                 :: Tv1            ! First Tangent Vector at CP
    DOUBLE PRECISION, DIMENSION (3)                 :: Tv2            ! Second Tangent Vector at CP
    ! Product is defined as [Tv1 x Tv2]

    JacM = matmul(Xa,Derived_PSF)    ! Calculates JacM from coordinates array and Derived_PSF
    Tv1 = JacM(:,1)                  ! Extracts Tv1
    Tv2 = JacM(:,2)                  ! Extracts Tv2

    !NormalVersor Calculation

    NormalVersor=cross_product(Tv1,Tv2)/norm2(cross_product(Tv1,Tv2))

    END SUBROUTINE CP_NormalVersor

    SUBROUTINE CP_Velocity (CP, FrstrmSpeed, OwnV)
    ! Stores the OwnVelocity and Wind Velocity at the Control Point.
    ! VELOCITIES DETERMINATION PROCEDURE IS MISSING!!!!!---------------

    !Input Variables
    DOUBLE PRECISION, dimension (3) ::    OwnV ! control point own velocity
    DOUBLE PRECISION, dimension (3) ::    FrstrmSpeed ! free stream velocity at control point

    !Output Variales
    type(Control_Point), INTENT (INOUT) :: CP

    CP%OwnVelocity = OwnV
    CP%WindVelocity = FrstrmSpeed


    END SUBROUTINE CP_Velocity

    END MODULE classControl_Point