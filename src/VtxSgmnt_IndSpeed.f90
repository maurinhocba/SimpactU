
    FUNCTION VtxSgmnt_IndSpeed (NP1, NP2, CP, Gamma)
    !Calculates the speed induced by a vortex segment (from NP1 to NP2) of intensity Gamma, in the control point CP.

    IMPLICIT NONE


    interface
    FUNCTION cross_product(U,V)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
    DOUBLE PRECISION, DIMENSION(3)              :: cross_product
    DOUBLE PRECISION                            :: cross1, cross2, cross3

    END FUNCTION cross_product
    end interface

    !Input Variables
    DOUBLE PRECISION, DIMENSION (3)  :: NP1, NP2    ! Nodal Points of the Segment, first-last.
    DOUBLE PRECISION, DIMENSION (3)  :: CP          ! Control Point to calculate induced speed in.
    DOUBLE PRECISION                 :: Gamma       ! Circulation of the segment.

    !Output Variables
    DOUBLE PRECISION, DIMENSION (3)  :: VtxSgmnt_IndSpeed   ! Induced speed in the control point.

    !Local Variables
    DOUBLE PRECISION, DIMENSION (3)  :: R1, R2, RL, R3        ! Pointing vectors from CP to NPi.
    DOUBLE PRECISION, PARAMETER      :: PI=3.1415926535897933 ! Pi.
    DOUBLE PRECISION                 :: L               ! Length of the vortex segment.
    DOUBLE PRECISION                 :: delta = 0.01    ! "Cutoff" radius parameter (VAN GARREL)
    DOUBLE PRECISION, DIMENSION (3)  :: e1, e2
    DOUBLE PRECISION                 :: NR1,NR2
    !Pointing Vectors
    R1 = CP-NP1
    R2 = CP-NP2
	
	RL = R1-R2
	R3 = cross_product(RL,R1) 

    !Norms
    NR1=norm2(R1)
    NR2=norm2(R2)

    ! Van Garrel Implementation of cutoff radius

    L=norm2(NP1-NP2) ! Length of the vortex segment

    !Phillips & Snyder Formula: Determines the induced speed at the Cp, uses Van Garrel cut off

    !VtxSgmnt_IndSpeed = Gamma*(NR1+NR2)*(cross_product(R1,R2))/(4.d0*PI*(NR1*NR2*(NR1*NR2+dot_product(R1,R2))+(delta*L)**2.d0))  !Parenthesis added to factor out (4pi)

 !   e1 = R1/norm2(R1)
 !   e2 = R2/norm2(R2)
 !   !
 !   if ((norm2(cross_product((R1-R2),R1))/(norm2(R1-R2))) >= delta*L) then
 !
 !       VtxSgmnt_IndSpeed = Gamma*(cross_product((R1-R2),R1))*dot_product((R1-R2),(e1-e2))/(4*PI*(norm2(cross_product((R1-R2),R1)))**2)
 !
 !   else
 !       VtxSgmnt_IndSpeed= [0.d0,0.d0,0.d0]
	!end if
		
    e1 = R1/NR1
    e2 = R2/NR2
    
    if ((norm2(R3)/(norm2(RL))) >= delta*L) then
    
        VtxSgmnt_IndSpeed = Gamma*R3*dot_product(RL,(e1-e2))/(4.d0*PI*(norm2(R3))**2.d0)
      
    else
        VtxSgmnt_IndSpeed= [0.d0,0.d0,0.d0]
    end if

    !end if

    !   if (VtxSgmnt_IndSpeed(1) < 1e-16) then
    !           VtxSgmnt_IndSpeed(1) = 0.0
    !       end if
    !
    !       if (VtxSgmnt_IndSpeed(2) < 1e-16) then
    !           VtxSgmnt_IndSpeed(2) = 0.0
    !       end if
    !
    !       if (VtxSgmnt_IndSpeed(3) < 1e-16) then
    !           VtxSgmnt_IndSpeed(3) = 0.0
    !       end if
    !if (VtxSgmnt_IndSpeed(1) /= VtxSgmnt_IndSpeed(1)) then
    !     print*,'NP', NP1, NP2
    !     print*,'VTX IS >> ', VtxSgmnt_IndSpeed
    !     print*,'L >> ', L
    !     pause
    !end if

    END FUNCTION VtxSgmnt_IndSpeed