
#include "fintrf.h"
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
    use celestialFunctions
    use JB2008Module
    implicit none

    mwPointer :: plhs(*), prhs(*)
    integer :: nlhs, nrhs

    mwPointer :: mxCreateDoubleScalar
    real(kind = 8) :: mxGetScalar
    real(kind = 8) :: julianDay, alt, lat, lon, f10, f10A, s10
    real(kind = 8) :: s10A, xm10, xm10A, y10, y10A, dstdtc, density
    real(kind = 8) :: sat(3), sun(2), temp(2)
    real(kind = 8) :: pi, degreesToRadians
    real(kind = 8) :: days1950, modifiedJulianDay
    real(kind = 8) :: sunRA, sunDec, gwRA

    pi = 4.0 * atan(1.0)
    degreesToRadians = pi / 180.0

    if(nrhs /= 13) then
       call mexErrMsgTxt('JB2008: 13 inputs required!')
    endif


    julianDay = mxGetScalar(prhs(1))
    alt = mxGetScalar(prhs(2))
    lat = mxGetScalar(prhs(3))
    lon = mxGetScalar(prhs(4))
    f10 = mxGetScalar(prhs(5))
    f10A = mxGetScalar(prhs(6))
    s10 = mxGetScalar(prhs(7))
    s10A = mxGetScalar(prhs(8))
    xm10 = mxGetScalar(prhs(9))
    xm10A = mxGetScalar(prhs(10))
    y10 = mxGetScalar(prhs(11))
    y10A = mxGetScalar(prhs(12))
    dstdtc = mxGetScalar(prhs(13))

    modifiedJulianDay = julianDay - 2400000.5
    days1950 = modifiedJulianDay - 33281.0

    gwRA = THETA(days1950)
    sat(1) = mod(gwRA + lon * degreesToRadians + 2.0 * pi, 2.0 * pi)
    sat(2) = lat * degreesToRadians
    sat(3) = alt

    call SUNPOS(modifiedJulianDay, sunRA, sunDec)
    sun(1) = sunRA
    sun(2) = sunDec

    call JB2008(modifiedJulianDay,sun,sat,f10,f10A,s10,s10A,xm10,xm10A, &
                y10,y10A,dstdtc,temp,density)

    plhs(1) = mxCreateDoubleScalar(temp(1))
    plhs(2) = mxCreateDoubleScalar(temp(2))
    plhs(3) = mxCreateDoubleScalar(density)

end subroutine mexfunction
