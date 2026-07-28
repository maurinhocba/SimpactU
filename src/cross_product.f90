    FUNCTION cross_product(U,V)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)  :: U,V
    DOUBLE PRECISION, DIMENSION(3)              :: cross_product
    DOUBLE PRECISION                            :: cross1, cross2, cross3

    cross1=U(2)*V(3)-U(3)*V(2)
    cross2=U(3)*V(1)-U(1)*V(3)
    cross3=U(1)*V(2)-U(2)*V(1)

    cross_product = [cross1,cross2,cross3]

    END FUNCTION cross_product


