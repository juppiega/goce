!>@file
!> DTM2013 Example program

program dtm2013_example_program

!This program is an example on how to use the DTM 2013 wrapper

use dtm_wrapper_interfaces  !Load the interface for the "dtm_wrapper" subroutine
use dtm2013_interfaces      !Load the interface for the "dtm2013", "P_ReadDTM12", and
                            !"density uncertainty" subroutines
use t_dtm_date              !Provides the dtm_date custom type

implicit none

!Define pi
real pi
parameter(pi=3.14159265)

!INPUT for dtm_wrapper
    type(dtm_date) today !the "dtm_date" type is defined within the t_dtm_date
                         !it provides a convenient way to input the universal time
                         !to the dtm wrapper routine

    real alti !altitude in km
    real alat !latitude in degrees
    real xlon !longitude in degrees

!OUTPUT from dtm_wrapper
      real :: ro     !density at the given altitude
      real :: tinf   !exospheric temperature
      real :: tz     !temperature at the given altitude
      real :: tp120  !vertical temperature gradient at 120 km
      real :: wmm    !mean molecular mass
      real :: ro_unc !Uncertainty in density
      real :: d(6)   !Concentrations of:
                     ! 1 atomic hydrogen
                     ! 2 helium
                     ! 3 atomic oxygen
                     ! 4 molecular nitrogen
                     ! 5 molecular oxygen
                     ! 6 atomic nitrogen (currently unused)

      !optional arguments for dtm_wrapper
      real, dimension(2) :: f_out
      real, dimension(2) :: fbar_out
      real, dimension(4) :: akp_out
      real :: hl_out
      real :: dayofyear_out

!Begin program

call load_config(784,"config.cfg")

!Set inputs
 alti=200.0 !real
 xlon=-4.0
 alat=40.0

! set a date

today%type_flag=2 ! Set this flag to 2 to input a calendar date
today%day=3
today%month=9
today%year=2011
today%hour=8
today%minute=32
today%second=34.23d0 ! warning! real*8


!call dtm wrapper

call dtm_wrapper(today,alti,alat,xlon,tz,tinf,tp120,ro,ro_unc,d,wmm, &
                 f_out=f_out,fbar_out=fbar_out,akp_out=akp_out, &
                 hl_out=hl_out,dayofyear_out=dayofyear_out)


!output the results

      write(*,*) " "
      write(*,'("-----------CALL TO DTM2013-----------")')
      write(*,*) " "


      write(*,*) "inputs"
      write(*,'("   Date:             ", f12.7, "  ( ",i2,"/",i2,"/",i4, " ) at ", i2, ":", i2, ":", f12.9)') today%mjd2000, &
                                         today%day, today%month,today%year, today%hour, today%minute, today%second
      write(*,'("   Day of year:      ", f10.6)') dayofyear_out
      write(*,'("   local time (hours): ", f12.7)') hl_out
      write(*,'("   latitude:         ", f12.7)') alat
      write(*,'("   longitude:        ", f12.7)') xlon
      write(*,'("   altitude:         ", f12.7)') alti
      write(*,'("   f :               ", 2(f12.7,1x))') f_out
      write(*,'("   fbar :            ", 2(f12.7,1x))') fbar_out
      write(*,'("   akp :             ", 4(f12.7,1x))') akp_out
      write(*,*) " "
      write(*,*) "outputs"
      write(*,'("   Temp at altitude :", f12.7)') tz
      write(*,'("   exospheric tmp :  ", f12.7)') tinf
      write(*,'("   atomic hydrogen : ", ES16.7)') d(1)
      write(*,'("   helium :          ", ES16.7)') d(2)
      write(*,'("   atomic oxygen :   ", ES16.7)') d(3)
      write(*,'("   molecular nitro : ", ES16.7)') d(4)
      write(*,'("   molecular oxygen :", ES16.7)') d(5)
      !write(*,'("   atomic nitrogen : ", ES16.7)') d(6) Ignore this value
      write(*,'("   density (g/cm^3) :", ES16.7, "  (+- ", f6.2, "% )")') ro, ro_unc
      write(*,'("   mean molec mass : ", ES16.7)') wmm

end program
