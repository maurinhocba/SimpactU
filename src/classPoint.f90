    MODULE classPoint

    ! "Point" abstract data type definition

    ! Martín Eduardo Pérez Segura
    ! Mauro S. Maza

    IMPLICIT NONE
    PRIVATE

    TYPE, Public :: Point
        !Attributes:
        DOUBLE PRECISION, dimension (3) ::   xyz ! point coordinates
        DOUBLE PRECISION, dimension (3) ::   xyz_loc !local coordinates

    CONTAINS ! methods' names

    END TYPE Point

    CONTAINS ! ==================================================================


    END MODULE classPoint