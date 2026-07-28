    MODULE SolverAXB

    !Solves the system [A] * {X} = {B}

    PUBLIC

    CONTAINS

    SUBROUTINE SOLVER (A_Mat, G, RHS)

    !Solves the system [A] * {G} = RHS, for the aerodynamic grid.

    !Input Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: A_Mat
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(IN)   :: RHS

    !Output Variables
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:), INTENT(OUT)   :: G

    !Local Variables
    INTEGER :: nPanels   !Size of the system, equals the number of grid panels.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: a  ! Dummy variable for function calling
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: b  ! Dummy variable for function calling
    integer, allocatable, dimension(:):: indx(:)   ! ludcmp variables
    DOUBLE PRECISION                              :: d               ! ludcmp variables
    integer                           :: np              ! ludcmp variables

    INTEGER, ALLOCATABLE, DIMENSION (:) :: IPIV
    INTEGER :: INFO

    nPanels = size(RHS) !Variable definition
    np = nPanels
    a = A_Mat           !Variable definition
    b = RHS             !Variable definition

    allocate(indx(nPanels))
    allocate (IPIV(np))

    call DGESV( np, 1, a, np, IPIV, b, np, INFO )

    ! Unknown vector replacement
    G = b;

    !deallocate(a, b)
    deallocate(IPIV)

    END SUBROUTINE SOLVER

    END MODULE SolverAXB
