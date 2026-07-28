    MODULE ShFunct_DShFunct

    ! Module containing Shape Functions and its derivatives.

    IMPLICIT NONE
    PUBLIC

    CONTAINS  !Module Procedures

    SUBROUTINE Shape_Functions (PanelType, PShapeFunct, XIv)
    !Returns the Shape Functions vector of a panel evaluated in XIv=(xi,etta,citta)

    !Input Variable
    INTEGER, INTENT(IN)            :: PanelType  ! Defines panel configuration according to internal code.
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: XIv ! Sets the computational point where the shape functions will be evaluated.

    !Output Variable
    DOUBLE PRECISION, DIMENSION (8), INTENT(OUT) :: PShapeFunct  ! Shape function vector of the panel evaluated in (XI,ETTA)=(0,0)

    !Local Variables
    DOUBLE PRECISION :: N1, N2, N3, N4, N5, N6, N7, N8
    DOUBLE PRECISION :: etta, xi, citta

    xi = XIv(1)
    etta = XIv(2)
    citta = XIv(3)

    N1= (etta-1)*(xi-1)*(citta+1)/4
    N2= -(etta-1)*(xi+1)*(citta+1)/4
    N3= (etta+1)*(xi+1)*(citta+1)/4
    N4= -(etta+1)*(xi-1)*(citta+1)/4
    N5= 1.d0/2
    N6= 1.d0/2
    N7= 1.d0/2
    N8= 1.d0/2

    select case (PanelType)

    case (3)
        PShapeFunct= [(1+xi+etta)*(citta+1)/3, (1-xi)*(citta+1)/3, (1+etta)*(citta+1)/3, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0]

    case (1)
        N5= 0.d0; N6= 0.d0; N7= 0.d0; N8= 0.d0
    case (51)
        N6= 0.d0; N7= 0.d0; N8= 0.d0
    case (52)
        N7= 0.d0; N8= 0.d0; N5= 0.d0
    case (53)
        N8= 0.d0; N5= 0.d0; N6= 0.d0
    case (54)
        N5= 0.d0; N6= 0.d0; N7= 0.d0
    case (61)
        N7= 0.d0; N8= 0.d0
    case (62)
        N8= 0.d0; N5= 0.d0
    case (63)
        N5= 0.d0; N6= 0.d0
    case (64)
        N6= 0.d0; N7= 0.d0
    case (71)
        N8= 0.d0
    case (72)
        N5= 0.d0
    case (73)
        N6= 0.d0
    case (74)
        N7= 0.d0
    end select

    if (PanelType /= 3) then
        PShapeFunct= [N1-(N5+N8)/2,N2-(N5+N6)/2,N3-(N6+N7)/2,N4-(N7+N8)/2,N5,N6,N7,N8]
    end if

    !print *,'psf : ',PShapeFunct
    END SUBROUTINE

    SUBROUTINE Derived_Shape_Functions (PanelType, Derived_PSF)
    !Returns the Shape Functions Derivatives vector of a panel evaluated in (XI,ETTA)=(0,0)

    !Input Variable
    INTEGER, INTENT(IN)   :: PanelType  ! Defines panel configuration according to internal code.

    !Output Variable
    DOUBLE PRECISION, DIMENSION (8,2), INTENT(OUT) :: Derived_PSF  ! Derived Shape function vector of the panel evaluated in (XI,ETTA)=(0,0)

    !Local Variables
    DOUBLE PRECISION :: N1e, N2e, N3e, N4e, N5e, N6e, N7e, N8e, N1s, N2s, N3s, N4s, N5s, N6s, N7s, N8s

    N1e= -1.d0/4.d0; N1s= -1.d0/4; N2e= 1.d0/4; N2s= -1.d0/4; N3e= 1.d0/4; N3s= 1.d0/4; N4e= -1.d0/4; N4s= 1.d0/4.d0
    N5e= 0.d0;    N5s= -1.d0/2; N6e= 1.d0/2; N6s= 0.d0;    N7e= 0.d0;   N7s= 1.d0/2; N8e= -1.d0/2; N8s= 0.d0

    select case (PanelType)

    case (3)

        Derived_PSF(:,1)= [1.d0/3, -1.d0/3, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
        Derived_PSF(:,2)= [1.d0/3, 0.d0, -1.d0/3, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0]

    case (1)
        N5e= 0.d0 ; N5s=0.d0;  N6e= 0.d0 ; N6s= 0.d0;  N7e= 0.d0 ; N7s= 0.d0;  N8e= 0.d0 ; N8s= 0.d0
    case (51)
        N6e= 0.d0 ; N6s= 0.d0;  N7e= 0.d0 ; N7s= 0.d0;  N8e= 0.d0 ; N8s= 0.d0
    case (52)
        N7e= 0.d0 ; N7s= 0.d0;  N8e= 0.d0 ; N8s= 0.d0;  N5e= 0.d0 ; N5s= 0.d0
    case (53)
        N8e= 0.d0 ; N8s= 0.d0;  N5e= 0.d0 ; N5s= 0.d0;  N6e= 0.d0 ; N6s= 0.d0
    case (54)
        N5e= 0.d0 ; N5s= 0.d0;  N6e= 0.d0 ; N6s= 0.d0;  N7e= 0.d0 ; N7s= 0.d0
    case (61)
        N7e= 0.d0 ; N7s= 0.d0;  N8e= 0.d0 ; N8s= 0.d0
    case (62)
        N8e= 0.d0 ; N8s= 0.d0;  N5e= 0.d0 ; N5s= 0.d0
    case (63)
        N5e= 0 ; N5s= 0.d0;  N6e= 0.d0 ; N6s= 0.d0
    case (64)
        N6e= 0.d0 ; N6s= 0.d0;  N7e= 0.d0 ; N7s= 0.d0
    case (71)
        N8e= 0.d0 ; N8s= 0.d0
    case (72)
        N5e= 0.d0 ; N5s= 0.d0
    case (73)
        N6e= 0.d0 ; N6s= 0.d0
    case (74)
        N7e= 0.d0 ; N7s= 0.d0
    end select

    if (PanelType /= 3) then
        Derived_PSF(:,1)= [N1e-(N5e+N8e)/2,N2e-(N5e+N6e)/2,N3e-(N6e+N7e)/2,N4e-(N7e+N8e)/2,N5e,N6e,N7e,N8e]
        Derived_PSF(:,2)= [N1s-(N5s+N8s)/2,N2s-(N5s+N6s)/2,N3s-(N6s+N7s)/2,N4s-(N7s+N8s)/2,N5s,N6s,N7s,N8s]
    end if

    END SUBROUTINE

    END MODULE ShFunct_DShFunct
