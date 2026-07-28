MODULE P_S_N

    !Module for input data processing.

    USE classNodal_Point
    !USE classPanel


    IMPLICIT NONE
    PUBLIC

    CONTAINS  !Module Procedures

SUBROUTINE Panel_Type (IENArray, PanelType)
    ! Generates an array with the internal codification of panel types of each panel

    !Input Variables
    INTEGER, DIMENSION (:,:), INTENT(IN) :: IENArray

    !Output Variables
    INTEGER, DIMENSION (:), ALLOCATABLE, INTENT(OUT) :: PanelType

    !Local Variables
    INTEGER                :: NumberOfPanels
    INTEGER                :: PnlType
    INTEGER, DIMENSION (8) :: P
    INTEGER                :: j !Loop indexes

    NumberOfPanels = size(IENArray,2)
    PnlType = 0

    ALLOCATE (PanelType (NumberOfPanels))

    !PanelType = 1  ! Trial version only ---------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!บบบบบบบบบบ

    do j=1,NumberOfPanels
        P = IENArray(:,j)
        If (P(4)==0) then
            PnlType = 3  ! 3-nodes Triangle
        else
            if (P(5)==0 .and. P(6)==0 .and. P(7)==0 .and. P(8)==0) then
                PnlType = 1 ! 4-nodes Cuad
            else
                if (P(6)==0 .and. P(7)==0 .and. P(8)==0) then
                    PnlType = 51 ! 5-nodes Cuad: extra node in posic 5
                end if
                if (P(5)==0 .and. P(7)==0 .and. P(8)==0) then
                    PnlType = 52 ! 5-nodes Cuad: extra node in posic 6
                end if
                if (P(5)==0 .and. P(6)==0 .and. P(8)==0) then
                    PnlType = 53 ! 5-nodes Cuad: extra node in posic 7
                end if
                if (P(5)==0 .and. P(6)==0 .and. P(7)==0) then
                    PnlType = 54 ! 5-nodes Cuad: extra node in posic 8
                end if
            end if
        end if

        if (PnlType == 0) then
            if (P(7)==0 .and. P(8)==0) then
                PnlType = 61 ! 6-nodes Cuad: extra node in posic 5,6
            end if
            if (P(8)==0 .and. P(5)==0) then
                PnlType = 62 ! 6-nodes Cuad: extra node in posic 6,7
            end if
            if (P(5)==0 .and. P(6)==0) then
                PnlType = 63 ! 6-nodes Cuad: extra node in posic 7,8
            end if
            if (P(6)==0 .and. P(7)==0) then
                PnlType = 64 ! 6-nodes Cuad: extra node in posic 8,5
            end if
        end if

        if (PnlType == 0) then
            if (P(8)==0) then
                PnlType = 71 ! 7-nodes Cuad: extra node in posic 5,6,7
            else if (P(5)==0) then
                PnlType = 72 ! 7-nodes Cuad: extra node in posic 6,7,8
            else if (P(6)==0) then
                PnlType = 73 ! 7-nodes Cuad: extra node in posic 7,8,5
            else if (P(7)==0) then
                PnlType = 74 ! 7-nodes Cuad: extra node in posic 8,5,6
            else
                PnlType = 8 ! 8-nodes Cuad
            end if

        end if

        PanelType(j) = PnlType

    end do

END SUBROUTINE Panel_Type

SUBROUTINE Segments (GridLabel, IENArray, PanelType, NumberOfSegments, ISNArray, IESArray,ShPanels_Links, ShPanels_Storage)
    ! Generates VtxSegmentSet for a grid, from IENArray and PanelType.

    !INPUT VARIABLES
    INTEGER                              :: GridLabel
    INTEGER, DIMENSION (:,:), INTENT(IN) :: IENArray
    INTEGER, DIMENSION (:), INTENT(IN)   :: PanelType

    !OUTPUT VARIABLES
    INTEGER, INTENT(OUT)                               :: NumberOfSegments
    INTEGER, ALLOCATABLE, DIMENSION (:,:), INTENT(OUT) :: ISNArray, IESArray
    INTEGER, ALLOCATABLE, DIMENSION (:), INTENT(OUT)   :: ShPanels_Links, ShPanels_Storage

    !LOCAL VARIABLES
    INTEGER                               :: SgmntIndex, ctrl, AUXindex, PnlType ! Internal Control Variables
    INTEGER                               :: NumberOfPanels
    INTEGER, ALLOCATABLE, DIMENSION (:,:) :: SgmntArrayTEMP,  SgmntArrayAUX ! Aux. arrays
    INTEGER, ALLOCATABLE, DIMENSION (:)   :: ShPanels_AUX ! Aux. arrays
    INTEGER, DIMENSION (4)                :: rot  !Aux. array

    INTEGER :: a,b,f,g,h,i,j,k,l,m,n,s,t !Loops indexes

    NumberOfPanels=size(IENArray,2)

    print *, 'NumberOfPanels: ',NumberOfPanels

    allocate (SgmntArrayAUX (3, 8*NumberOfPanels))
    allocate (SgmntArrayTEMP (2, 8*NumberOfPanels))
    allocate (IESArray (8,NumberOfPanels))
    SgmntArrayAUX = 0  !Initialize with zeros
    SgmntArrayTEMP = 0 !Initialize with zeros

    IESArray = 0   !Initialize with zeros
    AUXindex = 1   !Initialize in one

    !SgmntArrayAUX:
    do n=1,NumberOfPanels
        PnlType = PanelType(n)

        select case (PnlType)

        case (1)       ! 4-node Cuad
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (2,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (4,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (4,n)
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1

        case (3)        ! 3-node Triag
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (2,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (3,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1

        case (51, 52, 53, 54)       ! 5-node Cuad
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
            if (PnlType==51) then
                SgmntArrayAUX(2,AUXindex) = IENArray (5,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (2,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            if (PnlType==52) then
                SgmntArrayAUX(2,AUXindex) = IENArray (6,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (4,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            if (PnlType==53) then
                SgmntArrayAUX(2,AUXindex) = IENArray (7,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (4,n)
            if (PnlType==54) then
                SgmntArrayAUX(2,AUXindex) = IENArray (8,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            if (PnlType==51) then
                SgmntArrayAUX(1,AUXindex) = IENArray (5,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
            end if
            if (PnlType==52) then
                SgmntArrayAUX(1,AUXindex) = IENArray (6,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            end if
            if (PnlType==53) then
                SgmntArrayAUX(1,AUXindex) = IENArray (7,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            end if
            if (PnlType==54) then
                SgmntArrayAUX(1,AUXindex) = IENArray (8,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (4,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1

        case (61, 62, 63, 64)       !6-node Cuad
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
            if (PnlType==61 .or. PnlType==64) then
                SgmntArrayAUX(2,AUXindex) = IENArray (5,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (2,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            if (PnlType==61 .or. PnlType==62) then
                SgmntArrayAUX(2,AUXindex) = IENArray (6,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (4,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            if (PnlType==62) then
                SgmntArrayAUX(2,AUXindex) = IENArray (7,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (4,n)
            if (PnlType==63 .or. PnlType==64) then
                SgmntArrayAUX(2,AUXindex) = IENArray (8,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            if (PnlType==61 .or. PnlType==64) then
                SgmntArrayAUX(1,AUXindex) = IENArray (5,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if
            if (PnlType==61 .or. PnlType==62) then
                SgmntArrayAUX(1,AUXindex) = IENArray (6,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if
            if (PnlType==62 .or. PnlType==63) then
                SgmntArrayAUX(1,AUXindex) = IENArray (7,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if
            if (PnlType==63 .or. PnlType==64) then
                SgmntArrayAUX(1,AUXindex) = IENArray (8,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (4,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if

        case (71, 72, 73, 74, 8)       !7-node Cuad / 8-node Cuad
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
            if (PnlType==71 .or. PnlType==73 .or. PnlType==74 .or. PnlType==8) then
                SgmntArrayAUX(2,AUXindex) = IENArray (5,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (2,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            if (PnlType==71 .or. PnlType==72 .or. PnlType==74 .or. PnlType==8) then
                SgmntArrayAUX(2,AUXindex) = IENArray (6,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (4,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
            if (PnlType==71 .or. PnlType==72 .or. PnlType==73 .or. PnlType==8) then
                SgmntArrayAUX(2,AUXindex) = IENArray (7,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            SgmntArrayAUX(1,AUXindex) = IENArray (1,n)
            SgmntArrayAUX(2,AUXindex) = IENArray (4,n)
            if (PnlType==72 .or. PnlType==73 .or. PnlType==74 .or. PnlType==8) then
                SgmntArrayAUX(2,AUXindex) = IENArray (8,n)
            end if
            SgmntArrayAUX(3,AUXindex) = n
            AUXindex = AUXindex + 1
            if (PnlType==71 .or. PnlType==73 .or. PnlType==74 .or. PnlType==8) then
                SgmntArrayAUX(1,AUXindex) = IENArray (5,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (2,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if
            if (PnlType==71 .or. PnlType==72 .or. PnlType==74 .or. PnlType==8) then
                SgmntArrayAUX(1,AUXindex) = IENArray (6,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if
            if (PnlType==71 .or. PnlType==72 .or. PnlType==73 .or. PnlType==8) then
                SgmntArrayAUX(1,AUXindex) = IENArray (7,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (3,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if
            if (PnlType==72 .or. PnlType==73 .or. PnlType==74 .or. PnlType==8) then
                SgmntArrayAUX(1,AUXindex) = IENArray (8,n)
                SgmntArrayAUX(2,AUXindex) = IENArray (4,n)
                SgmntArrayAUX(3,AUXindex) = n
                AUXindex = AUXindex + 1
            end if

        end select
    end do
    print*, 'AUXindex: ', AUXindex, '(/)'
    !print*, SgmntArrayAUX
    !------
    !SgmntArrayTEMP:
    !IES:

    SgmntIndex = 0

    do n=1, AUXindex  ! Trial version only: was  NumberOfPanels*8 !!!!!!!!--------บบบบบบ
        do m=1,(n-1)
            if (((SgmntArrayAUX(1,n) == SgmntArrayTEMP(1,m) .and. SgmntArrayAUX(2,n) == SgmntArrayTEMP(2,m))) .or. &
                ((SgmntArrayAUX(1,n) == SgmntArrayTEMP(2,m) .and. SgmntArrayAUX(2,n) == SgmntArrayTEMP(1,m)))) then
                do i=1,8
                    if (IESArray (i,SgmntArrayAUX(3,n)) == 0) then
                        IESArray(i,SgmntArrayAUX(3,n)) = m
                        exit
                    end if
                end do
                exit
            end if
        end do

        if ((SgmntArrayAUX(1,n) /= SgmntArrayTEMP(1,m) .or. SgmntArrayAUX(2,n) /= SgmntArrayTEMP(2,m)) .and. &
            (SgmntArrayAUX(1,n) /= SgmntArrayTEMP(2,m) .and. SgmntArrayAUX(2,n) /= SgmntArrayTEMP(1,m))) then

            SgmntIndex = SgmntIndex + 1
            SgmntArrayTEMP(1:2,SgmntIndex) = SgmntArrayAUX(1:2,n)
            do i=1,8

                if (IESArray (i,SgmntArrayAUX(3,n)) == 0) then
                    IESArray(i,SgmntArrayAUX(3,n)) = SgmntIndex
                    exit
                end if
                !exit
            end do
        end if
	end do

    NumberOfSegments = SgmntIndex

    !IES segment local number relocation

    do j=1,NumberOfPanels
        rot = 0 !Initialize with zeros
        if (PanelType(j) == 72) then
            rot(2:4) = IESArray(5:7,j)
        end if
        if (PanelType(j) == 73) then
            rot(3:4) = IESArray(5:6,j)
            rot(1) = IESArray(7,j)
        end if
        if (PanelType(j) == 74) then
            rot(1:2) = IESArray(6:7,j)
            rot(4) = IESArray(5,j)
        end if
        if (PanelType(j) == 62) then
            rot(2:3) = IESArray(5:6,j)
        end if
        if (PanelType(j) == 63) then
            rot(3:4) = IESArray(5:6,j)
        end if
        if (PanelType(j) == 64) then
            rot(1) = IESArray(6,j)
            rot(4) = IESArray(5,j)
        end if
        if (PanelType(j) == 52) then
            rot(2) = IESArray(5,j)
        end if
        if (PanelType(j) == 53) then
            rot(3) = IESArray(5,j)
        end if
        if (PanelType(j) == 54) then
            rot(4) = IESArray(5,j)
        end if

        IESArray(5:8,j) = rot(:)
    end do

    print *, 'NumberOfSegments_PSN: ', NumberOfSegments
    !------
    !ISNArray:

    allocate (ISNArray (2, NumberOfSegments))

    ISNArray(:,:) = SgmntArrayTEMP(:,1:NumberOfSegments)

    ! print *, 'ISNArray: '
    !do n=1,NumberOfSegments
    !    print *, n, ' : ',ISNArray(:,n),'/'
    !end do
    !pause
    !print *, 'IESArray: '
    !do n=1,NumberOfPanels
    !    print *, n, ' : ',IESArray(:,n),'/'
    !end do
    !pause

    !------
    !Sharing Panels
    allocate (ShPanels_AUX (NumberOfSegments))
    allocate (ShPanels_Links (NumberOfSegments))
    allocate (ShPanels_Storage (4*NumberOfPanels))
    ShPanels_AUX = 0      !Initialize with zeros
    ShPanels_Links = 0    !Initialize with zeros
    ShPanels_Storage = 0  !Initialize with zeros

    do f=1,NumberOfPanels
        do g=1,8
            if (IESArray(g,f)/= 0) then
                ShPanels_AUX(IESArray(g,f))=ShPanels_AUX(IESArray(g,f)) + 1
            end if
        end do
    end do

    do h=2,NumberOfSegments
        ShPanels_AUX(h)=ShPanels_AUX(h)+ShPanels_AUX(h-1)
    end do

    ShPanels_Links=ShPanels_AUX ! Links vector
    !print *, 'ShPanels_Links: ', ShPanels_Links,'//' !!!!!!!!!!!!!!!!!!!!!!!!! PRINTING CHECK

    do a=1,NumberOfPanels
        do b=1,8
            if (IESArray(b,a)/= 0) then
                ShPanels_Storage(ShPanels_AUX(IESArray(b,a)))=a
                ShPanels_AUX(IESArray(b,a))=ShPanels_AUX(IESArray(b,a))-1
            end if
        end do
    end do

    !print *, 'ShPanels_Storage: ', ShPanels_Storage
    !pause

    !Memory deallocation
    deallocate (SgmntArrayTEMP)
    deallocate(SgmntArrayAUX)
    deallocate(ShPanels_AUX)


    END SUBROUTINE Segments

END MODULE P_S_N
