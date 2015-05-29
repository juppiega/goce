! *****************************************************************
! MATLAB® mex-interface for the HWM07 model.                      *
!                                                                 *
! By Juho Iipponen (Finnish Meteorological Institute)             *
!                                                                 *
! Building:                                                       *
! mex FCFLAGS="\$FCFLAGS -O3" -output hwm07_mex apexcord.F90 dwm07b.F90 HWM07Module.F90 hwm07_mex.F90                                                      *
!                                                                 *
! Usage:                                                          *
! [zonalW, meridionalW]                                           *
! = hwm07_mex(year, day_of_year, altitude, latitude, longitude,ap)*
! *****************************************************************

#include "fintrf.h"
!#include "mex.h"
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

    use HWM07Module ! The model functions.

    implicit none

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
    real :: apNow  ! Geomagnetic ap proxy
    real :: year, doy      ! Year and day of year (0-366, decimal)
    real :: ap(2) ! Model input array. ap(1) = 0, ap(2) = apNow
    character(128) :: defaultdata = 'hwm071308e.dat' ! Input file name.
    integer :: yearDayInput, twoDigitYear
    real :: secondOfDay

    ! Outputs
    real :: wind(2) ! Model result. wind(1) = Meridional, wind(2) = Zonal

    real :: pi, zero = 0.0

    character(len = 80) temp

    pi = 4.0 * atan(1.0)

    
    ! Check input and output argument counts
    if(nrhs /= 6 .and. nrhs /= 0) then
       call mexErrMsgTxt('HWM07: 6 or 0 inputs required!')
    endif

    if(nrhs == 0) then
        call loadmodel(defaultdata)
        
        zero = 0.0
        plhs(1) = mxCreateDoubleScalar(dble(zero))
        plhs(2) = mxCreateDoubleScalar(dble(zero))

        return
    endif

    if(nlhs > 2) then
       call mexErrMsgTxt('HWM07: Too many outputs!')
    endif

    ! Read input values
    year = mxGetScalar(prhs(1))
    doy = mxGetScalar(prhs(2))
    alt = mxGetScalar(prhs(3))
    lat = mxGetScalar(prhs(4))
    lon = mxGetScalar(prhs(5))
    apNow = mxGetScalar(prhs(6))

    ! Convert year and doy into correct input types.
    twoDigitYear = year - 100 * floor(year / 100.0)
    yearDayInput = 1000 * twoDigitYear + floor(doy)
    secondOfDay = (doy - floor(doy)) * 86400.0

    ! Convert the scalars into correct input arrays
    ap(1) = 0
    ap(2) = apNow

    ! Print Debug information
    !write(temp, *) f(1)
    !k=mexPrintf('f(1): '//temp//achar(13))
    !write(temp, *) fbar(1)
    !k=mexPrintf('fbar(1): '//temp//achar(13))
    !write(temp, *) km(1)
    !k=mexPrintf('km(1): '//temp//achar(13))
    !write(temp, *) yearDayInput
    !k=mexPrintf('iyd: '//temp//achar(13))
    !write(temp, *) secondOfDay
    !k=mexPrintf('sec: '//temp//achar(13))
    !write(temp, *) alt
    !k=mexPrintf('alt: '//temp//achar(13))
    !write(temp, *) lat
    !k=mexPrintf('lat: '//temp//achar(13))
    !write(temp, *) lon
    !k=mexPrintf('lon: '//temp//achar(13))
    !write(temp, *) ap(2)
    !k=mexPrintf('ap(2): '//temp//achar(13))

    ! Now call HWM07 **************************************************************
    call hwm07(yearDayInput, secondOfDay, alt, lat, lon, zero,zero,zero,ap,wind) !*
    ! *****************************************************************************

    !write(temp, *) wind(2)
    !k=mexPrintf('u: '//temp//achar(13))

    plhs(1) = mxCreateDoubleScalar(dble(wind(2)))
    plhs(2) = mxCreateDoubleScalar(dble(wind(1)))
    
    return
end subroutine mexfunction
