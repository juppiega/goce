! *****************************************************************
! MATLAB® mex-interface for the DTM-2013 model.                   *
!                                                                 *
! Adapted by Juho Iipponen (Finnish Meteorological Institute)     *
! from an example program "dtm2013_example_direct_call"           *
! available in the DTM-2013 user manual pdf.                      *
!                                                                 *
! Building:                                                       *
! mex FCFLAGS="\$FCFLAGS -O3" -output dtm2013_mex t_dtm_date.F90 dtm_interfaces.F90 DTM_2012_subroutines.F90 dtm2013_mex.F90                                                             *
!                                                                 *
! Usage:                                                          *
! [tinf, tz, rho, rho_unc, wmm, d]                                *
! = dtm2013_mex(day_of_year, altitude, latitude, longitude,       *
!            local_solar_time, F10, F10A, Am3hAgo, Am24hAverage)  *
! *****************************************************************

#include "fintrf.h"
!#include "mex.h"
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

    use dtm2013_interfaces     ! Load the interface for the "dtm2013", "P_ReadDTM13", and
                               ! "density uncertainty" subroutines

    implicit none

    interface a2K
        function a2K(a_input)
            real, intent(in) :: a_input 
        end function
    end interface

    mwPointer :: plhs(*), prhs(*) ! Matlab input and output arguments
    integer :: nlhs, nrhs

    ! Matlab utility functions for converting between Fortran and Matlab formats
    mwPointer :: mxCreateDoubleScalar
    mwPointer :: mxGetPr
    mwPointer :: mxCreateDoubleMatrix
    mwPointer :: output_ptr
    mwSize    :: output_size
    real(kind = 8) :: mxGetScalar
    integer(kind = 4) :: mexPrintf
    character(len = 80) :: mess
    integer :: k

    ! Inputs
    real :: alt, lat, lon ! Position (km and degrees)
    real :: f10, f10A     ! Solar proxies
    real :: km3h, km24hA, am3h, am24hA  ! Geomagnetic am proxies
    real :: lst, doy      ! Local solar time (hours, decimal) and day of year (0-366, decimal)
    real :: f(2), fbar(2), km(4) ! Model input arrays

    ! Outputs
    real :: density, tinf, tz, tg120 ! Density, exosph.- and local temperature,
                                     ! Temperature gradient at 120 km.
    real :: meanMolarMass, densityUncert
    real :: d(6) ! Concentrations
                 ! 1 atomic hydrogen
                 ! 2 helium
                 ! 3 atomic oxygen
                 ! 4 molecular nitrogen
                 ! 5 molecular oxygen
                 ! 6 atomic nitrogen (currently unused)

    real :: pi, degreesToRadians, m, c

    pi = 4.0 * atan(1.0)
    degreesToRadians = pi / 180.0
    m = 0.521660985377791 ! Constants to convert am -> km
    c = 1.575168678069221
    
    ! Check input and output argument counts
    if(nrhs /= 9) then
       call mexErrMsgTxt('DTM2013: 9 inputs required!')
    endif

    if(nlhs > 6) then
       call mexErrMsgTxt('DTM2013: Too many outputs!')
    endif

    ! Read input values
    doy = mxGetScalar(prhs(1))
    alt = mxGetScalar(prhs(2))
    lat = mxGetScalar(prhs(3)) * degreesToRadians
    lon = mxGetScalar(prhs(4)) * degreesToRadians
    lst = (mxGetScalar(prhs(5)) / 24.0) * 2.0 * pi
    f10 = mxGetScalar(prhs(6))
    f10A = mxGetScalar(prhs(7))
    am3h = mxGetScalar(prhs(8))
    am24hA = mxGetScalar(prhs(9))

    !Convert am -> km
    km3h = max(0.0, a2K(am3h))
    km24hA = max(0.0, a2K(am24hA))

    ! Convert the scalars into correct input arrays
    f(1) = f10
    f(2) = 0.0
    fbar(1) = f10A
    fbar(2) = 0.0
    km(1) = km3h
    km(2) = 0.0
    km(3) = km24hA
    km(4) = 0.0

    !First, a call to P_ReadDTM13 is required to initialize the dtm2013 model
    call P_ReadDTM12 ()

    !Now call dtm2013
    call dtm2012(doy, f, fbar, km, alt, lst, lat, lon, tz, tinf, tg120, &
                 density, d, meanMolarMass)
    
    ! Convert lst back to hours, and compute uncertainty
    !lst = lst * 24.0 / (2.0 * pi)
    !call density_uncertainty(alt, lat, lst, f(1), km(1), densityUncert)
    
    densityUncert = 1E-30

    plhs(1) = mxCreateDoubleScalar(dble(tinf))
    plhs(2) = mxCreateDoubleScalar(dble(tz))
    plhs(3) = mxCreateDoubleScalar(dble(density))
    plhs(4) = mxCreateDoubleScalar(dble(densityUncert))
    plhs(5) = mxCreateDoubleScalar(dble(meanMolarMass))
    
    ! The sixth output is a matrix of concentrations: d
    plhs(6) = mxCreateDoubleMatrix(6,1,0) ! 6 by 1 real matrix
    output_ptr = mxGetPr(plhs(6))
    output_size = 6
    call mxCopyReal8ToPtr(dble(d), output_ptr, output_size) 

    return
end subroutine mexfunction

real function a2K(a_input)
    implicit none

    !input variable
    real, intent(in) :: a_input

    !Conversion tables
    integer :: table_length
    parameter (table_length=28)


    real :: AM(table_length), AKM(table_length)

    DATA AM/  0.  ,  1.4 ,  3.4 ,  5.4 ,  7.4 , 10.4 , 13.4 , 16.4 , &
             20.4 , 26.4 , 33.4 , 40.4 , 50.4 , 60.4 , 70.4 , 86.4 , &
            103.4 ,120.4 ,146.4 ,173.4 ,200.4 ,243.4 ,286.4 ,330.4 , &
            386.4 ,443.4 ,500.4 ,611.4/

    DATA AKM/ 0.  ,  0.33,  0.66,  1.  ,  1.33,  1.66,  2.  ,  2.33, &
              2.66,  3.  ,  3.33,  3.66,  4.  ,  4.33,  4.66,  5.  , &
              5.33,  5.66,  6.  ,  6.33,  6.66,  7.  ,  7.33,  7.66, &
              8.  ,  8.33,  8.66,  9./

    !Auxiliaries
    integer :: i, im
    real :: a


    !a might be modified inside the routine.
	!Define an internal a, so we can keep the intent(in) for a_input
    a=a_input


    a2K = -1.0
    do  i=1,table_length

              if(a.EQ.AM(i)) then
                  a2K=AKM(i)
                  exit
              else if (a.LE.AM(i)) then
                  im=i-1
                  a2K=AKM(im) + (AKM(i) - AKM(im))*(a - AM(im)) / (AM(i) - AM(im))
                  exit
              end if
     end do


end function
