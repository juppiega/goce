#include "fintrf.h"

! Fit IL-Model using the Levenberg-Marquardt method.
! Building:
! mex -O FCFLAGS="\$FCFLAGS -std=f2008" -output levenbergMarquardt_mex lmSolver.F90 levenbergMarquardt_mex.F90 -llapack

module fitModule
    implicit none
    
    type dataStruct
        real(kind = 8), allocatable, dimension(:) :: P10, P11, P20, P21, P22, &
                                     P30, P31, P32, P33, P40, P41, P42, P43, P44, P50, P51, &
                                     P52, P53, P54, P60, P62, P63, P64, P70, P71, P74, &
                                     mP10,mP20,mP30,mP40,mP50,mP60,mP70,mP11,mP31,mP51, &
                                     mP21, mP41, mP32, mP52, yv, dv, latitudeTerm, solarTerm, &
                                     annual, lv, dv_mag, &
                                     diurnal, semidiurnal, terdiurnal, quaterdiurnal,&
                                     geomagnetic, data, Z, F, FA
        integer(kind = 4), allocatable :: coeffInd(:)
        integer :: name, numBiases
        real(kind = 8) :: dT0, T0
        real(kind = 8), allocatable :: aeInt(:,:), biases(:,:)
    end type

    type(dataStruct) :: TexStruct, OStruct, N2Struct, HeStruct, ArStruct, O2Struct, rhoStruct
    real(kind = 8), allocatable :: weights(:), initGuess(:), dTCoeffs(:), T0Coeffs(:)

    interface clamp
        module procedure clamp_scalar, clamp_vector
    end interface

contains

function structToDerived_TexAndMajor(matlabStruct, typeName) result(D)
    implicit none
    ! Input
    mwPointer, intent(in) :: matlabStruct
    character(len = *), intent(in) :: typeName
    ! Output
    type(dataStruct) :: D
    ! Function declarations
    mwPointer mxGetField
    mwPointer mxGetM, mxGetN
    ! Local variables
    mwPointer mrows, ncols, field, mxGetPr
    mwSize N, M, numel
    real(kind = 8) :: numBiases_real
    real(kind = 8), allocatable :: coeffInd_real(:), aeInt_temp(:)
    character(len = 80) tempChar
    integer(kind = 4) :: mexPrintf, k

    ! Find out size
    mrows = mxGetM(mxGetField(matlabStruct, 1, 'data'))
    N = mrows

    ! Allocate
    allocate(D%P10(N), D%P11(N), D%mP10(N), D%mP11(N))
    allocate(D%P20(N), D%P21(N), D%P22(N), D%mP20(N), D%mP21(N))
    allocate(D%P30(N), D%P31(N), D%P32(N), D%P33(N), D%mP30(N), D%mP31(N), D%mP32(N))
    allocate(D%P40(N), D%P41(N), D%P42(N), D%P43(N), D%P44(N), D%mP40(N), D%mP41(N))
    allocate(D%P50(N), D%P51(N), D%P52(N), D%P53(N), D%P54(N),D%mP50(N),D%mP51(N),D%mP52(N))
    allocate(D%P60(N), D%P62(N), D%P63(N), D%P64(N), D%mP60(N))
    allocate(D%P70(N), D%P71(N), D%P74(N), D%mP70(N))
    allocate(D%yv(N), D%dv(N), D%lv(N), D%dv_mag(N), D%latitudeTerm(N), D%solarTerm(N))
    allocate(D%annual(N), D%diurnal(N), D%semidiurnal(N), D%terdiurnal(N))
    allocate(D%quaterdiurnal(N), D%geomagnetic(N), D%data(N), D%Z(N), D%F(N), D%FA(N))

    !write(tempChar,*) size(D%data)
    !k = mexPrintf('size(data): '//tempChar//achar(13))
    
    
    ! Copy information
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P10')), D%P10, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P11')), D%P11, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP10')), D%mP10, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP11')), D%mP11, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P20')), D%P20, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP20')), D%mP20, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP21')), D%mP21, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P21')), D%P21, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P22')), D%P22, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P30')), D%P30, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P31')), D%P31, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP30')), D%mP30, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP31')), D%mP31, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP32')), D%mP32, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P32')), D%P32, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P33')), D%P33, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P40')), D%P40, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P41')), D%P41, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P42')), D%P42, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP40')), D%mP40, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP41')), D%mP41, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P43')), D%P43, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P44')), D%P44, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P50')), D%P50, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P51')), D%P51, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP50')), D%mP50, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP51')), D%mP51, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP52')), D%mP52, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P52')), D%P52, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P53')), D%P53, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P54')), D%P54, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P60')), D%P60, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP60')), D%mP60, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P62')), D%P62, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P63')), D%P63, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P64')), D%P64, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P70')), D%P70, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'mP70')), D%mP70, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P71')), D%P71, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'P74')), D%P74, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'yv')), D%yv, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'dv')), D%dv, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'dv_mag')), D%dv_mag, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'lv')), D%lv, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'latitudeTerm')), D%latitudeTerm, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'solarTerm')), D%solarTerm, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'annual')), D%annual, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'diurnal')), D%diurnal, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'semidiurnal')), D%semidiurnal, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'terdiurnal')), D%terdiurnal, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'quaterdiurnal')), D%quaterdiurnal, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'geomagnetic')), D%geomagnetic, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'data')), D%data, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'Z')), D%Z, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'F')), D%F, N)
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'FA')), D%FA, N)

    ! AE Integral
    mrows = mxGetM(mxGetField(matlabStruct, 1, 'aeInt'))
    ncols = mxGetN(mxGetField(matlabStruct, 1, 'aeInt'))
    M = mrows; N = ncols
    numel = mrows * ncols
    write(tempChar, *) numel
    allocate(D%aeInt(M,N))
    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'aeInt')), D%aeInt, numel)
    
    if (typeName /= 'rho') then
        ! Coefficient indices
        mrows = max(mxGetM(mxGetField(matlabStruct, 1, 'coeffInd')), mxGetN(mxGetField(matlabStruct, 1, 'coeffInd')))
        N = mrows
        allocate(D%coeffInd(N), coeffInd_real(N))
        call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'coeffInd')), coeffInd_real, N)
        D%coeffInd = nint(coeffInd_real)
    end if

    mrows = mxGetM(mxGetField(matlabStruct, 1, 'data'))
    N = mrows
    if (typeName == 'O') then
        D%name = 1
    else if (typeName == 'N2') then
        D%name = 2
    else if (typeName == 'He') then
        D%name = 3
    else if (typeName == 'Ar') then
        D%name = 4
    else if (typeName == 'O2') then
        D%name = 5
        D%numBiases = 0
    end if
    if (typename == 'O' .or. typename == 'N2' .or. typename == 'He' .or. &
	typename == 'Ar') then
	    call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'numBiases')), numBiases_real, 1)
        D%numBiases = nint(numBiases_real)
        allocate(D%biases(N, D%numBiases))
        call mxCopyPtrToReal8(mxGetPr(mxGetField(matlabStruct, 1, 'biases')), D%biases, D%numBiases*N)
    endif
    !write(tempChar, *) D%numBiases
    !k = mexPrintf('numBiases '//tempChar//achar(13))

    return
end function

subroutine deallocateStruct(D, typename)
    implicit none
    type(dataStruct), intent(inout) :: D
    character(len = *), intent(in) :: typename

    deallocate(D%P10, D%P11, D%mP10, D%mP11)
    deallocate(D%P20, D%P21, D%P22, D%mP20, D%mP21)
    deallocate(D%P30, D%P31, D%P32, D%P33, D%mP30, D%mP31, D%mP32)
    deallocate(D%P40, D%P41, D%P42, D%P43, D%P44, D%mP40, D%mP41)
    deallocate(D%P50, D%P51, D%P52, D%P53, D%P54, D%mP50, D%mP51, D%mP52)
    deallocate(D%P60, D%P62, D%P63, D%P64, D%mP60)
    deallocate(D%P70, D%P71, D%P74, D%mP70)
    deallocate(D%yv, D%dv, D%lv, D%dv_mag, D%latitudeTerm, D%solarTerm)
    deallocate(D%annual, D%diurnal, D%semidiurnal, D%terdiurnal)
    deallocate(D%quaterdiurnal, D%geomagnetic, D%data, D%Z, D%F, D%FA)
    deallocate(D%aeInt)

    if (typename /= 'rho') then
        deallocate(D%coeffInd)
    end if
    if (typename == 'O' .or. typename == 'N2' .or. typename == 'He' .or. &
	typename == 'Ar') then
	    deallocate(D%biases)
	end if
end subroutine

function G_major(a, S, numBiases)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    integer, intent(in) :: numBiases
    real(kind = 8), allocatable :: G_major(:), latitudeTerm(:), solarTerm(:), annual(:), diurnal(:), semidiurnal(:), &
                                   terdiurnal(:), quaterdiurnal(:), geomagnetic(:), geom_symmetric(:), geom_yearly(:),&
                                   geom_lst(:), AE_base(:), geom_solar(:), geom_lon(:), longitudinal(:)
    integer :: k, dPh, numInts
    real(kind=8) :: dPy
    integer(kind = 4) :: mexPrintf, i
    real(kind = 8) :: pi
    pi = 4.0 * atan(1.0)
    
    k = numBiases + 1; ! Counter, which helps adding termS%
    !i = mexPrintf('G Begin'//achar(13))
    ! Latitude termS%
    latitudeTerm = a(k+1)*S%P10 + a(k+2)*S%P20 + a(k+3)*S%P30 + a(k+4)*S%P40 + a(k+5)*S%P50 + a(k+6)*S%P60 + &
                     a(k+7)*S%FA*S%P10 + a(k+8)*S%FA*S%P20 + a(k+9)*S%FA*S%P30 + a(k+10)*S%FA*S%P40 + &
                     a(k+11)*S%F*S%P10 + a(k+12)*S%F*S%P20 + a(k+13)*S%F*S%P30 + a(k+14)*S%F*S%P40;
    k = k + 14;

    ! % Solar activity termS%
    !i = mexPrintf('G:solar'//achar(13))
    solarTerm = a(k+1)*S%F + a(k+2)*S%F**2 + a(k+3)*S%FA + a(k+4)*S%FA**2 + a(k+5)*S%F*S%FA;
    k = k + 5;

    ! Annual termS%
    !i = mexPrintf('G:annual'//achar(13))
    annual = (a(k+1) + a(k+2)*S%P10 + a(k+3)*S%P20 + a(k+4)*S%P30 + a(k+5)*S%P40 + a(k+6)*S%FA + a(k+7)*S%F) * &
               (a(k+8)*sin(S%yv) + a(k+9)*cos(S%yv) + a(k+10)*sin(2*S%yv) + a(k+11)*cos(2*S%yv) + a(k+12)*sin(3*S%yv) + &
                a(k+13)*cos(3*S%yv) + a(k+14)*sin(4*S%yv) + a(k+15)*cos(4*S%yv));
    k = k + 15;

    dPh = k + 11;
    !i = mexPrintf('G:diurnal'//achar(13))
    diurnal = ((a(k+1)*S%P11 + a(k+2)*S%P31 + a(k+3)*S%P51 + a(k+4)*S%P71 + a(k+5)*S%FA + a(k+6)*S%FA**2 + &
                 a(k+7)*(S%F - S%FA)) + (a(k+8)*S%P11 + a(k+9)*S%P21 + a(k+10)*S%P31)*(cos(S%yv-pi*a(dPh))))*cos(S%dv) + &
                ((a(k+12)*S%P11 + a(k+13)*S%P31 + a(k+14)*S%P51 + a(k+15)*S%P71 + a(k+16)*S%FA + a(k+17)*S%FA**2 + &
                a(k+18)*(S%F - S%FA)) + (a(k+19)*S%P11 + a(k+20)*S%P21 + a(k+21)*S%P31)*(cos(S%yv-pi*a(dPh))))*sin(S%dv);
    k = k + 21;
    
    !i = mexPrintf('G:semidiurnal'//achar(13))
    semidiurnal = (a(k+1)*S%P22 + a(k+2)*S%P32 + a(k+3)*S%P52 + a(k+4)*S%FA + a(k+5)*S%FA**2 + a(k+6)*(S%F-S%FA) + &
                (a(k+7)*S%P32 + a(k+8)*S%P52)*cos(S%yv-pi*a(dPh)))*cos(2*S%dv) + &
                (a(k+9)*S%P22 + a(k+10)*S%P32 + a(k+11)*S%P52 + a(k+12)*S%FA + a(k+13)*S%FA**2 + a(k+14)*(S%F-S%FA) + &
                (a(k+15)*S%P32 + a(k+16)*S%P52)*cos(S%yv-pi*a(dPh)))*sin(2*S%dv);
    k = k + 16;

    !i = mexPrintf('G:terdiurnal'//achar(13))
    terdiurnal = (a(k+1)*S%P33 + a(k+2)*S%P53 + (a(k+3)*S%P43 + a(k+4)*S%P63)*cos(S%yv-pi*a(dPh)))*cos(3*S%dv) + &
                   (a(k+5)*S%P33 + a(k+6)*S%P53 + (a(k+7)*S%P43 + a(k+8)*S%P63)*cos(S%yv-pi*a(dPh)))*sin(3*S%dv);
    k = k + 8;

    !i = mexPrintf('G:quaterdiurnal'//achar(13))
    quaterdiurnal = a(k+1)*S%P44*cos(4*S%dv) + a(k+2)*S%P44*sin(4*S%dv);
    k = k + 2;

    longitudinal = (1.0 + a(k+1)*S%P44 + a(k+2)*S%P64)*cos(S%dv - 4*S%lv + a(k+3))*fourier4(S, a(k+4:k+8))+&
                 (1.0 + a(k+9)*S%P44 + a(k+10)*S%P64)*cos(2*S%dv - 4*S%lv + a(k+11))*fourier4(S, a(k+12:k+16))+&
                 (1.0 + a(k+17)*S%P33 + a(k+18)*S%P53)*cos(S%dv - 3*S%lv + a(k+19))*fourier4(S, a(k+20:k+24))+&
                 (1.0 + a(k+25)*S%P22 + a(k+26)*S%P42)*cos(2*S%dv - 2*S%lv + a(k+27))*fourier4(S, a(k+28:k+32))+&
                 (1.0 + a(k+33)*S%FA)*(a(k+34)*S%P21+a(k+35)*S%P31+a(k+36)*S%P41)*cos(S%lv)+&
                 (1.0 + a(k+37)*S%FA)*(a(k+38)*S%P21+a(k+39)*S%P31+a(k+40)*S%P41)*sin(S%lv);
    k = k + 40;

    numInts = size(S%aeInt, 2)
    !i = mexPrintf('G:geomagnetic'//achar(13))
    
    !i = mexPrintf('G:AE_base'//achar(13))
    AE_base = sumRowWise([a(k+1), a(k+2), a(k+3), a(k+4), a(k+5), a(k+6), a(k+7)], S%aeInt)
    geom_symmetric = (a(k+8) + a(k+9)*S%mP20 + a(k+10)*S%mP40 + a(k+11)*S%mP60)*AE_base;
    dPy = a(k+16);
    geom_lon = (a(k+12) + a(k+13)*S%P21 + a(k+14)*S%P41)*(1+a(k+15)*cos(S%yv-pi*dPy))*AE_base*cos(S%lv-pi*a(k+17)) + &
           (a(k+18) + a(k+19)*S%P32 + a(k+20)*S%P52)*(1+a(k+21)*cos(S%yv-pi*dPy))*AE_base*cos(2*(S%lv-pi*a(k+22)));
    geom_lst = (a(k+23) + a(k+24)*S%P21 + a(k+25)*S%P41)*(1+a(k+26)*cos(S%yv-pi*dPy))*AE_base*cos(S%dv-pi*a(k+27)) + &
           (a(k+28) + a(k+29)*S%P32 + a(k+30)*S%P52)*(1+a(k+31)*cos(S%yv-pi*dPy))*AE_base*cos(2*(S%dv-pi*a(k+32)));
    geom_solar = (a(k+33) + a(k+34)*S%P10*cos(S%yv-pi*dPy) + a(k+35)*S%mP20)*AE_base*S%FA;
    geomagnetic = geom_symmetric + geom_lon + geom_lst + geom_solar;

    k = k + 21;

    !i = mexPrintf('G:sum'//achar(13))
    G_major = latitudeTerm + solarTerm + annual + diurnal + semidiurnal + terdiurnal + quaterdiurnal + geomagnetic;
    !i = mexPrintf('G End'//achar(13))
end function

function G_Tex(a, S, numBiases)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    integer, intent(in) :: numBiases
    real(kind = 8), allocatable :: G_Tex(:), latitudeTerm(:), solarTerm(:), annual(:), diurnal(:), semidiurnal(:), &
                                   terdiurnal(:), quaterdiurnal(:), geomagnetic(:), geom_symmetric(:), geom_yearly(:),&
                                   geom_lst(:), AE_base(:), geom_solar(:), geom_lon(:), longitudinal(:)
    integer :: k, dPh, numInts
    real(kind=8) :: dPy
    integer(kind = 4) :: mexPrintf, i
    real(kind = 8) :: pi
    pi = 4.0 * atan(1.0)
    
    k = numBiases + 1; ! Counter, which helps adding termS%
    !i = mexPrintf('G Begin'//achar(13))
    ! Latitude termS%
    latitudeTerm = a(k+1)*S%P10 + a(k+2)*S%P20 + a(k+3)*S%P30 + a(k+4)*S%P40 + a(k+5)*S%P50 + a(k+6)*S%P60 + &
                     a(k+7)*S%FA*S%P10 + a(k+8)*S%FA*S%P20 + a(k+9)*S%FA*S%P30 + a(k+10)*S%FA*S%P40 + &
                     a(k+11)*S%F*S%P10 + a(k+12)*S%F*S%P20 + a(k+13)*S%F*S%P30 + a(k+14)*S%F*S%P40;
    k = k + 14;

    ! % Solar activity termS%
    !i = mexPrintf('G:solar'//achar(13))
    solarTerm = a(k+1)*S%F + a(k+2)*S%F**2 + a(k+3)*S%FA + a(k+4)*S%FA**2 + a(k+5)*S%F*S%FA;
    k = k + 5;

    ! Annual termS%
    !i = mexPrintf('G:annual'//achar(13))
    annual = (a(k+1) + a(k+2)*S%P10 + a(k+3)*S%P20 + a(k+4)*S%P30 + a(k+5)*S%P40 + a(k+6)*S%FA + a(k+7)*S%F) * &
               (a(k+8)*sin(S%yv) + a(k+9)*cos(S%yv) + a(k+10)*sin(2*S%yv) + a(k+11)*cos(2*S%yv) + a(k+12)*sin(3*S%yv) + &
                a(k+13)*cos(3*S%yv) + a(k+14)*sin(4*S%yv) + a(k+15)*cos(4*S%yv));
    k = k + 15;

    dPh = k + 11;
    !i = mexPrintf('G:diurnal'//achar(13))
    diurnal = ((a(k+1)*S%P11 + a(k+2)*S%P31 + a(k+3)*S%P51 + a(k+4)*S%P71 + a(k+5)*S%FA + a(k+6)*S%FA**2 + &
                 a(k+7)*(S%F - S%FA)) + (a(k+8)*S%P11 + a(k+9)*S%P21 + a(k+10)*S%P31)*(cos(S%yv-pi*a(dPh))))*cos(S%dv) + &
                ((a(k+12)*S%P11 + a(k+13)*S%P31 + a(k+14)*S%P51 + a(k+15)*S%P71 + a(k+16)*S%FA + a(k+17)*S%FA**2 + &
                a(k+18)*(S%F - S%FA)) + (a(k+19)*S%P11 + a(k+20)*S%P21 + a(k+21)*S%P31)*(cos(S%yv-pi*a(dPh))))*sin(S%dv);
    k = k + 21;
    
    !i = mexPrintf('G:semidiurnal'//achar(13))
    semidiurnal = (a(k+1)*S%P22 + a(k+2)*S%P32 + a(k+3)*S%P52 + a(k+4)*S%FA + a(k+5)*S%FA**2 + a(k+6)*(S%F-S%FA) + &
                (a(k+7)*S%P32 + a(k+8)*S%P52)*cos(S%yv-pi*a(dPh)))*cos(2*S%dv) + &
                (a(k+9)*S%P22 + a(k+10)*S%P32 + a(k+11)*S%P52 + a(k+12)*S%FA + a(k+13)*S%FA**2 + a(k+14)*(S%F-S%FA) + &
                (a(k+15)*S%P32 + a(k+16)*S%P52)*cos(S%yv-pi*a(dPh)))*sin(2*S%dv);
    k = k + 16;

    !i = mexPrintf('G:terdiurnal'//achar(13))
    terdiurnal = (a(k+1)*S%P33 + a(k+2)*S%P53 + (a(k+3)*S%P43 + a(k+4)*S%P63)*cos(S%yv-pi*a(dPh)))*cos(3*S%dv) + &
                   (a(k+5)*S%P33 + a(k+6)*S%P53 + (a(k+7)*S%P43 + a(k+8)*S%P63)*cos(S%yv-pi*a(dPh)))*sin(3*S%dv);
    k = k + 8;

    !i = mexPrintf('G:quaterdiurnal'//achar(13))
    quaterdiurnal = a(k+1)*S%P44*cos(4*S%dv) + a(k+2)*S%P44*sin(4*S%dv);
    k = k + 2;

    dPh = k + 14;
    numInts = size(S%aeInt, 2)
    !i = mexPrintf('G:geomagnetic'//achar(13))
    
    longitudinal = (1.0 + a(k+1)*S%P44 + a(k+2)*S%P64)*cos(S%dv - 4*S%lv + a(k+3))*fourier4(S, a(k+4:k+8))+&
                 (1.0 + a(k+9)*S%P44 + a(k+10)*S%P64)*cos(2*S%dv - 4*S%lv + a(k+11))*fourier4(S, a(k+12:k+16))+&
                 (1.0 + a(k+17)*S%P33 + a(k+18)*S%P53)*cos(S%dv - 3*S%lv + a(k+19))*fourier4(S, a(k+20:k+24))+&
                 (1.0 + a(k+25)*S%P22 + a(k+26)*S%P42)*cos(2*S%dv - 2*S%lv + a(k+27))*fourier4(S, a(k+28:k+32))+&
                 (1.0 + a(k+33)*S%FA)*(a(k+34)*S%P21+a(k+35)*S%P31+a(k+36)*S%P41)*cos(S%lv)+&
                 (1.0 + a(k+37)*S%FA)*(a(k+38)*S%P21+a(k+39)*S%P31+a(k+40)*S%P41)*sin(S%lv);
    k = k + 40;

    numInts = size(S%aeInt, 2)
    !i = mexPrintf('G:geomagnetic'//achar(13))
    
    !i = mexPrintf('G:AE_base'//achar(13))
    AE_base = sumRowWise([a(k+1), a(k+2), a(k+3), a(k+4), a(k+5), a(k+6), a(k+7)], S%aeInt)
    geom_symmetric = (a(k+8) + a(k+9)*S%mP20 + a(k+10)*S%mP40 + a(k+11)*S%mP60)*AE_base;
    dPy = a(k+16);
    geom_lon = (a(k+12) + a(k+13)*S%P21 + a(k+14)*S%P41)*(1+a(k+15)*cos(S%yv-pi*dPy))*AE_base*cos(S%lv-pi*a(k+17)) + &
           (a(k+18) + a(k+19)*S%P32 + a(k+20)*S%P52)*(1+a(k+21)*cos(S%yv-pi*dPy))*AE_base*cos(2*(S%lv-pi*a(k+22)));
    geom_lst = (a(k+23) + a(k+24)*S%P21 + a(k+25)*S%P41)*(1+a(k+26)*cos(S%yv-pi*dPy))*AE_base*cos(S%dv-pi*a(k+27)) + &
           (a(k+28) + a(k+29)*S%P32 + a(k+30)*S%P52)*(1+a(k+31)*cos(S%yv-pi*dPy))*AE_base*cos(2*(S%dv-pi*a(k+32)));
    geom_solar = (a(k+33) + a(k+34)*S%P10*cos(S%yv-pi*dPy) + a(k+35)*S%mP20)*AE_base*S%FA;
    geomagnetic = geom_symmetric + geom_lon + geom_lst + geom_solar;

    k = k + 21;

    !i = mexPrintf('G:sum'//achar(13))
    G_Tex = latitudeTerm + solarTerm + annual + diurnal + semidiurnal + terdiurnal + quaterdiurnal + geomagnetic;
    !i = mexPrintf('G End'//achar(13))
end function

function G_lbT0(a, S, numBiases)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    integer, intent(in) :: numBiases
    real(kind = 8), allocatable :: G_lbT0(:), latitudeTerm(:), solarTerm(:), annual(:), diurnal(:), semidiurnal(:), &
                                   terdiurnal(:), quaterdiurnal(:), geomagnetic(:), geom_symmetric(:), geom_yearly(:),&
                                   geom_lst(:), AE_base(:)
    integer :: k, dPy, numInts
    integer(kind = 4) :: mexPrintf, i
    real(kind = 8) :: pi
    pi = 4.0 * atan(1.0)
    
    k = numBiases + 1;

    latitudeTerm = a(k+1)*S%P10 + a(k+2)*S%P20 + a(k+3)*S%P30 + a(k+4)*S%P40 + &
                     a(k+5)*S%FA*S%P10 + a(k+6)*S%FA*S%P20 + a(k+7)*S%FA*S%P30 + &
                     a(k+8)*S%F*S%P10 + a(k+9)*S%F*S%P20 + a(k+10)*S%F*S%P30;
    k = k + 10;

    solarTerm = a(k+1)*S%F + a(k+2)*S%F**2 + a(k+3)*S%FA + a(k+4)*S%FA**2 + a(k+5)*S%F*S%FA;
    k = k + 5;

    annual = (a(k+1) + a(k+2)*S%P10 + a(k+3)*S%P20 + a(k+4)*S%P30 + a(k+5)*S%P40 + a(k+6)*S%FA + a(k+7)*S%F) * &
               (a(k+8)*sin(S%yv) + a(k+9)*cos(S%yv) + a(k+10)*sin(2*S%yv) + a(k+11)*cos(2*S%yv));
    k = k + 11;

    dPy = k + 7;
    diurnal = ((a(k+1)*S%P11 + a(k+2)*S%P31 + a(k+3)*S%FA + a(k+4)*S%FA**2) + (a(k+5)*S%P11 + a(k+6)*S%P21)*&
                (cos(S%yv-pi*a(dPy))))*cos(S%dv) + &
                ((a(k+8)*S%P11 + a(k+9)*S%P31 + a(k+10)*S%FA + a(k+11)*S%FA**2) + (a(k+12)*S%P11 + a(k+13)*S%P21)*&
                (cos(S%yv-pi*a(dPy))))*sin(S%dv);
    k = k + 13;

    semidiurnal = ((a(k+1)*S%P22 + a(k+2)*S%P32 + a(k+3)*S%FA + a(k+4)*S%FA**2) + &
                    (a(k+5)*S%P32)*cos(S%yv-pi*a(dPy)))*cos(2*S%dv) + &
                    ((a(k+6)*S%P22 + a(k+7)*S%P32 + a(k+8)*S%FA + a(k+9)*S%FA**2) &
                    + (a(k+10)*S%P32)*cos(S%yv-pi*a(dPy)))*sin(2*S%dv);
    k = k + 10;

    !i = mexPrintf('G:sum'//achar(13))
    G_lbT0 = latitudeTerm + solarTerm + annual + diurnal + semidiurnal;
    !i = mexPrintf('G End'//achar(13))
end function

function G_minor(a, S, numBiases)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    integer, intent(in) :: numBiases
    real(kind = 8), allocatable :: G_minor(:), latitudeTerm(:), solarTerm(:), annual(:), diurnal(:), semidiurnal(:), &
                                   terdiurnal(:), quaterdiurnal(:), geomagnetic(:), geom_symmetric(:), geom_yearly(:),&
                                   geom_lst(:), AE_base(:)
    integer :: k, dPy, numInts
    integer(kind = 4) :: mexPrintf, i
    real(kind = 8) :: pi
    pi = 4.0 * atan(1.0)
    
    k = numBiases + 1;

    latitudeTerm = a(k+1)*S%P10 + a(k+2)*S%P20 + a(k+3)*S%P30 + a(k+4)*S%P40 + &
                     a(k+5)*S%FA*S%P10 + a(k+6)*S%FA*S%P20 + a(k+7)*S%FA*S%P30 + &
                     a(k+8)*S%F*S%P10 + a(k+9)*S%F*S%P20 + a(k+10)*S%F*S%P30;
    k = k + 10;

    solarTerm = a(k+1)*S%F + a(k+2)*S%F**2 + a(k+3)*S%FA + a(k+4)*S%FA**2 + a(k+5)*S%F*S%FA;
    k = k + 5;

    annual = (a(k+1) + a(k+2)*S%P10 + a(k+3)*S%P20 + a(k+4)*S%P30 + a(k+5)*S%P40 + a(k+6)*S%FA + a(k+7)*S%F) * &
               (a(k+8)*sin(S%yv) + a(k+9)*cos(S%yv) + a(k+10)*sin(2*S%yv) + a(k+11)*cos(2*S%yv));
    k = k + 11;

    dPy = k + 7;
    diurnal = ((a(k+1)*S%P11 + a(k+2)*S%P31 + a(k+3)*S%FA + a(k+4)*S%FA**2) + (a(k+5)*S%P11 + a(k+6)*S%P21)*&
                (cos(S%yv-pi*a(dPy))))*cos(S%dv) + &
                ((a(k+8)*S%P11 + a(k+9)*S%P31 + a(k+10)*S%FA + a(k+11)*S%FA**2) + (a(k+12)*S%P11 + a(k+13)*S%P21)*&
                (cos(S%yv-pi*a(dPy))))*sin(S%dv);
    k = k + 13;

    semidiurnal = ((a(k+1)*S%P22 + a(k+2)*S%P32 + a(k+3)*S%FA + a(k+4)*S%FA**2) + &
                    (a(k+5)*S%P32)*cos(S%yv-pi*a(dPy)))*cos(2*S%dv) + &
                    ((a(k+6)*S%P22 + a(k+7)*S%P32 + a(k+8)*S%FA + a(k+9)*S%FA**2) &
                    + (a(k+10)*S%P32)*cos(S%yv-pi*a(dPy)))*sin(2*S%dv);
    k = k + 10;

    !i = mexPrintf('G:sum'//achar(13))
    G_minor = latitudeTerm + solarTerm + annual + diurnal + semidiurnal;
    !i = mexPrintf('G End'//achar(13))
end function

function G_lbDT(a, S, numBiases)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    integer, intent(in) :: numBiases
    real(kind = 8), allocatable :: G_lbDT(:), latitudeTerm(:), solarTerm(:), annual(:), diurnal(:), semidiurnal(:), &
                                   terdiurnal(:), quaterdiurnal(:), geomagnetic(:), geom_symmetric(:), geom_yearly(:),&
                                   geom_lst(:), AE_base(:)
    integer :: k, dPh, numInts
    integer(kind = 4) :: mexPrintf, i
    real(kind = 8) :: pi
    pi = 4.0 * atan(1.0)
    
    k = numBiases + 1; ! Counter, which helps adding termS%
    !i = mexPrintf('G Begin'//achar(13))
    ! Latitude termS%
    latitudeTerm = a(k+1)*S%P10 + a(k+2)*S%P20 + a(k+3)*S%P30 + a(k+4)*S%P40 + a(k+5)*S%P50 + a(k+6)*S%P60 + &
                     a(k+7)*S%FA*S%P10 + a(k+8)*S%FA*S%P20 + a(k+9)*S%FA*S%P30 + a(k+10)*S%FA*S%P40 + &
                     a(k+11)*S%F*S%P10 + a(k+12)*S%F*S%P20 + a(k+13)*S%F*S%P30 + a(k+14)*S%F*S%P40;
    k = k + 14;

    ! % Solar activity termS%
    !i = mexPrintf('G:solar'//achar(13))
    solarTerm = a(k+1)*S%F + a(k+2)*S%F**2 + a(k+3)*S%FA + a(k+4)*S%FA**2 + a(k+5)*S%F*S%FA;
    k = k + 5;

    ! Annual termS%
    !i = mexPrintf('G:annual'//achar(13))
    annual = (a(k+1) + a(k+2)*S%P10 + a(k+3)*S%P20 + a(k+4)*S%P30 + a(k+5)*S%P40 + a(k+6)*S%FA + a(k+7)*S%F) * &
               (a(k+8)*sin(S%yv) + a(k+9)*cos(S%yv) + a(k+10)*sin(2*S%yv) + a(k+11)*cos(2*S%yv) + a(k+12)*sin(3*S%yv) + &
                a(k+13)*cos(3*S%yv) + a(k+14)*sin(4*S%yv) + a(k+15)*cos(4*S%yv));
    k = k + 15;

    dPh = k + 11;
    !i = mexPrintf('G:diurnal'//achar(13))
    diurnal = ((a(k+1)*S%P11 + a(k+2)*S%P31 + a(k+3)*S%P51 + a(k+4)*S%P71 + a(k+5)*S%FA + a(k+6)*S%FA**2 + &
                 a(k+7)*(S%F - S%FA)) + (a(k+8)*S%P11 + a(k+9)*S%P21 + a(k+10)*S%P31)*(cos(S%yv-pi*a(dPh))))*cos(S%dv) + &
                ((a(k+12)*S%P11 + a(k+13)*S%P31 + a(k+14)*S%P51 + a(k+15)*S%P71 + a(k+16)*S%FA + a(k+17)*S%FA**2 + &
                a(k+18)*(S%F - S%FA)) + (a(k+19)*S%P11 + a(k+20)*S%P21 + a(k+21)*S%P31)*(cos(S%yv-pi*a(dPh))))*sin(S%dv);
    k = k + 21;
    
    !i = mexPrintf('G:semidiurnal'//achar(13))
    semidiurnal = (a(k+1)*S%P22 + a(k+2)*S%P32 + a(k+3)*S%P52 + a(k+4)*S%FA + a(k+5)*S%FA**2 + a(k+6)*(S%F-S%FA) + &
                (a(k+7)*S%P32 + a(k+8)*S%P52)*cos(S%yv-pi*a(dPh)))*cos(2*S%dv) + &
                (a(k+9)*S%P22 + a(k+10)*S%P32 + a(k+11)*S%P52 + a(k+12)*S%FA + a(k+13)*S%FA**2 + a(k+14)*(S%F-S%FA) + &
                (a(k+15)*S%P32 + a(k+16)*S%P52)*cos(S%yv-pi*a(dPh)))*sin(2*S%dv);
    k = k + 16;

    !i = mexPrintf('G:terdiurnal'//achar(13))
    terdiurnal = (a(k+1)*S%P33 + a(k+2)*S%P53 + (a(k+3)*S%P43 + a(k+4)*S%P63)*cos(S%yv-pi*a(dPh)))*cos(3*S%dv) + &
                   (a(k+5)*S%P33 + a(k+6)*S%P53 + (a(k+7)*S%P43 + a(k+8)*S%P63)*cos(S%yv-pi*a(dPh)))*sin(3*S%dv);
    k = k + 8;

    !i = mexPrintf('G:quaterdiurnal'//achar(13))
    quaterdiurnal = a(k+1)*S%P44*cos(4*S%dv) + a(k+2)*S%P44*sin(4*S%dv);
    k = k + 2;

    !i = mexPrintf('G:sum'//achar(13))
    G_lbDT = latitudeTerm + solarTerm + annual + diurnal + semidiurnal + terdiurnal + quaterdiurnal
    !i = mexPrintf('G End'//achar(13))
end function

function fourier4(S, a)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    real(kind = 8) :: fourier4
    
    fourier4 = a(1) + a(2)*sin(S%yv) + a(3)*cos(S%yv) + a(4)*sin(2*S%yv) + a(5)*cos(2*S%yv)

end function

function evalTex(S, coeff)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: coeff(:)
    real(kind = 8), allocatable :: evalTex(:)
    integer(kind = 4) :: mexPrintf, k
    
    allocate(evalTex(size(S%data)))
    !k = mexPrintf('Before G at evalTex'//achar(13))
    evalTex = coeff(1) * (1 + clamp(dble(-0.9), G_Tex(coeff, S, 0), dble(9.0))) 

end function

function evalT0(S, coeff)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: coeff(:)
    real(kind = 8), allocatable :: evalT0(:)
    integer(kind = 4) :: mexPrintf, k
    
    allocate(evalT0(size(S%data)))
    !k = mexPrintf('Before G at evalT0'//achar(13))
    evalT0 = clamp(dble(300.0), coeff(1) * (1 + clamp(dble(-0.5), G_lbT0(coeff, S, 0), dble(2.0))), dble(1000.0))

end function

function evalDT(S, coeff)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: coeff(:)
    real(kind = 8), allocatable :: evalDT(:)
    integer(kind = 4) :: mexPrintf, k
    
    allocate(evalDT(size(S%data)))
    !k = mexPrintf('Before G at evalDT'//achar(13))
    evalDT = coeff(1) * (1 + clamp(dble(-0.9), G_lbDT(coeff, S, 0), dble(9.0))); 

end function

function evalMajorSpecies(S, coeff, numBiases)
    implicit none
    type(dataStruct), intent(in) :: S
    integer, intent(in) :: numBiases
    real(kind = 8), intent(in) :: coeff(:)
    real(kind = 8), allocatable :: evalMajorSpecies(:)
    integer(kind = 4) :: mexPrintf, k
    !character(len = 60) :: tempChar
    
    allocate(evalMajorSpecies(size(S%data)))
    !write(tempChar,*) S%numBiases
    !k = mexPrintf('S%numBiases: '//tempChar//achar(13))
    evalMajorSpecies = exp(coeff(1) + G_major(coeff, S, numBiases));

end function

function evalMinorSpecies(S, coeff, numBiases)
    implicit none
    type(dataStruct), intent(in) :: S
    integer, intent(in) :: numBiases
    real(kind = 8), intent(in) :: coeff(:)
    real(kind = 8), allocatable :: evalMinorSpecies(:)
    integer(kind = 4) :: mexPrintf, k
    !character(len = 60) :: tempChar
    
    allocate(evalMinorSpecies(size(S%data)))
    !write(tempChar,*) S%numBiases
    !k = mexPrintf('S%numBiases: '//tempChar//achar(13))
    evalMinorSpecies = exp(coeff(1) + G_minor(coeff, S, numBiases));

end function

subroutine computeDensityRHS(S, Tex, dT0, T0, rhs)
    implicit none
    type(dataStruct), intent(in) :: S
    real(kind = 8), intent(in) :: Tex(:), dT0(:), T0(:)
    real(kind = 8), allocatable, intent(out) :: rhs(:)
    integer :: numBiasesOrig
    real(kind = 8), parameter :: u2kg = 1.660538921E-27, k = 1.38064852E-23, g = 9.418
    real(kind = 8) :: alpha, molecMass
    real(kind = 8), dimension(size(Tex)) :: sigma, T, gamma, altTerm

    if (S%name == 1) then
        molecMass = 16 * u2kg;
        alpha = 0;
    elseif (S%name == 2) then
        molecMass = 28 * u2kg;
        alpha = 0;
    elseif (S%name == 3) then
        molecMass = 4 * u2kg;
        alpha = -0.38;
    elseif (S%name == 4) then
        molecMass = 40 * u2kg;
        alpha = 0;
    elseif (S%name == 5) then
        molecMass = 32 * u2kg;
        alpha = 0;
    end if

    sigma = dT0 / (Tex - T0);
    T = Tex - (Tex - T0) * exp(-sigma * (S%Z));
    gamma = molecMass * g / (k * sigma * 1E-3 * Tex);
    altTerm = (1 + gamma + alpha) * log(T0 / T) - gamma * sigma * (S%Z);
    rhs = log(S%data) - altTerm;

end subroutine

function computeRho(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens) result(rho)
    implicit none
    real(kind = 8), intent(in) :: Tex(:), dT0(:), T0(:), Z(:), olbDens(:), N2lbDens(:), HelbDens(:), &
                                  ArlbDens(:), O2lbDens(:)
    real(kind = 8), dimension(size(Tex)) :: sigma, T, gamma_O, gamma_N2, gamma_He, gamma_Ar, &
                                            gamma_O2,f_O,f_N2,f_He,f_Ar,f_O2,OnumDens,N2numDens, &
                                            HeNumDens, ArNumDens, O2NumDens, rho
    real(kind = 8), parameter :: u2kg = 1.660538921E-27, k = 1.38064852E-23, g = 9.418

    sigma = dT0 / (Tex - T0);

    T = Tex - (Tex - T0) * exp(-sigma * (Z));

    gamma_O = 16 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_O = (T0 / T)**(1+gamma_O) * exp(-sigma * (Z) * gamma_O);
    OnumDens = OlbDens*f_O; ! [1/cm^3]

    gamma_N2 = 28 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_N2 = (T0 / T)**(1+gamma_N2) * exp(-sigma * (Z) * gamma_N2);
    N2numDens = N2lbDens*f_N2; ! [1/cm^3]

    gamma_He = 4 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_He = (T0 / T)**(1+gamma_He-0.38) * exp(-sigma * (Z) * gamma_He);
    HeNumDens = HelbDens*f_He; ! [1/cm^3]
    
    gamma_Ar = 40 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_Ar = (T0 / T)**(1+gamma_Ar) * exp(-sigma * (Z) * gamma_Ar);
    ArNumDens = ArlbDens*f_Ar; ! [1/cm^3]

    gamma_O2 = 32 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_O2 = (T0 / T)**(1+gamma_O2) * exp(-sigma * (Z) * gamma_O2);
    O2NumDens = O2lbDens*f_O2; ! [1/cm^3]

    rho = (16*OnumDens + 28*N2numDens + 4*HeNumDens + 40*ArNumDens + 32*O2NumDens) * u2kg * 1E6; ! [kg/m^3]

end function

subroutine findTempsForFit(varStruct, TexStruct, coeff, dTCoeffs, T0Coeffs, Tex, dT0, T0)
    implicit none
    type(dataStruct), intent(in) :: varStruct, TexStruct
    real(kind = 8), intent(in) :: coeff(:), dTCoeffs(:), T0Coeffs(:)
    real(kind = 8), allocatable, intent(out) :: Tex(:), dT0(:), T0(:)
    real(kind = 8), allocatable :: Tex_est(:)

    Tex_est = evalTex(varStruct, coeff(TexStruct%coeffInd));

    T0 = evalT0(varStruct, T0Coeffs);
    dT0 = evalDT(varStruct, dTCoeffs);
    Tex = clamp(T0+1, Tex_est, dble(5000));

end subroutine

function computeMajorSpeciesResidual(varStruct, Tex, dT0, T0, coeff) result(residual)
    implicit none
    type(dataStruct) :: varStruct
    real(kind = 8), intent(in) :: Tex(:), dT0(:), T0(:), coeff(:)
    real(kind = 8), allocatable :: residual(:), Gvec(:), rhs(:)
    
    allocate(rhs(size(Tex)))
    call computeDensityRHS(varStruct, Tex, dT0, T0, rhs);
    Gvec = G_major(coeff, varStruct, varStruct%numBiases);

    if (varStruct%numBiases == 0) then
        residual = (rhs / max(coeff(1) + Gvec, dble(1))) - 1;
    elseif (varStruct%numBiases > 0) then
        residual = (rhs / max(coeff(1) + sumRowWise(coeff(2:varStruct%numBiases+1), varStruct%biases) + Gvec, dble(1)))&
                     - 1;
    end if

end function

function computeMinorSpeciesResidual(varStruct, Tex, dT0, T0, coeff) result(residual)
    implicit none
    type(dataStruct) :: varStruct
    real(kind = 8), intent(in) :: Tex(:), dT0(:), T0(:), coeff(:)
    real(kind = 8), allocatable :: residual(:), Gvec(:), rhs(:)
    
    allocate(rhs(size(Tex)))
    call computeDensityRHS(varStruct, Tex, dT0, T0, rhs);
    Gvec = G_minor(coeff, varStruct, varStruct%numBiases);

    if (varStruct%numBiases == 0) then
        residual = (rhs / max(coeff(1) + Gvec, dble(1))) - 1;
    elseif (varStruct%numBiases > 0) then
        residual = (rhs / max(coeff(1) + sumRowWise(coeff(2:varStruct%numBiases+1), varStruct%biases) + Gvec, dble(1)))&
                     - 1;
    end if

end function

function computeO2Residual(varStruct, Tex, dT0, T0, coeff) result(residual)
    implicit none
    type(dataStruct) :: varStruct
    real(kind = 8), intent(in) :: Tex(:), dT0(:), T0(:), coeff(:)
    real(kind = 8), allocatable :: residual(:), rhs(:)
    real(kind = 8) :: Gvec
    
    allocate(rhs(size(Tex)))
    call computeDensityRHS(varStruct, Tex, dT0, T0, rhs);
    Gvec = coeff(1);

    residual = (rhs / max(coeff(1) + Gvec, dble(1))) - 1;

end function

function modelMinimizationFunction(coeff) result(residual)
    implicit none
    real(kind = 8), intent(in) :: coeff(:)
    real(kind = 8), allocatable :: residual(:), TexMesuredEstimate(:), Tex(:), dT0(:),& 
                                   olbDens(:), N2lbDens(:), HelbDens(:), ArlbDens(:), O2lbDens(:),&
                                   modelRho(:), T0(:)
    integer(kind = 8), allocatable :: residInd(:)
    integer(kind = 8) :: dataLen, i, residIndEnd
    integer(kind = 4) :: mexPrintf, k, mexCallMatlab

    
    dataLen = size(TexStruct%data) + size(OStruct%data) + size(N2Struct%data) + size(HeStruct%data) + &
              size(rhoStruct%data) + size(ArStruct%data) + size(O2Struct%data)
    allocate(residual(dataLen))
    
    !k = mexPrintf('Before TexEst.'//achar(13))
    T0 = evalT0(TexStruct, T0Coeffs)
    TexMesuredEstimate = clamp(T0+1, evalTex(TexStruct, coeff(TexStruct%coeffInd)), dble(5000));
    !k = mexCallMatlab(0, 0, 1, TexMesuredEstimate(1), 'disp')
    residInd = (/(i, i = 1, size(TexStruct%data))/); residIndEnd = residInd(size(residInd))
    residual(residInd) = TexStruct%data/TexMesuredEstimate - 1;
    !deallocate(residInd)

    !k = mexPrintf('Before O'//achar(13))
    call findTempsForFit(OStruct, TexStruct, coeff, dTCoeffs, T0Coeffs, Tex, dT0, T0)
    residInd = residIndEnd + (/(i, i = 1, size(OStruct%data))/); residIndEnd = residInd(size(residInd))
    residual(residInd) = computeMajorSpeciesResidual(OStruct, Tex, dT0, T0, coeff(OStruct%coeffInd));
    !deallocate(residInd)

    !k = mexPrintf('Before N2'//achar(13))
    call findTempsForFit(N2Struct, TexStruct, coeff, dTCoeffs, T0Coeffs, Tex, dT0, T0);
    residInd = residIndEnd + (/(i, i = 1, size(N2Struct%data))/); residIndEnd = residInd(size(residInd))
    residual(residInd) = computeMajorSpeciesResidual(N2Struct, Tex, dT0, T0, coeff(N2Struct%coeffInd));
    !deallocate(residInd)

    !k = mexPrintf('Before He'//achar(13))
    call findTempsForFit(HeStruct, TexStruct, coeff, dTCoeffs, T0Coeffs, Tex, dT0, T0);
    residInd = residIndEnd + (/(i, i = 1, size(HeStruct%data))/); residIndEnd = residInd(size(residInd))
    residual(residInd) = computeMajorSpeciesResidual(HeStruct, Tex, dT0, T0, coeff(HeStruct%coeffInd));
    !deallocate(residInd)
    
    !k = mexPrintf('Before Ar'//achar(13))
    call findTempsForFit(ArStruct, TexStruct, coeff, dTCoeffs, T0Coeffs, Tex, dT0, T0);
    residInd = residIndEnd + (/(i, i = 1, size(ArStruct%data))/); residIndEnd = residInd(size(residInd))
    residual(residInd) = computeMinorSpeciesResidual(ArStruct, Tex, dT0, T0, coeff(ArStruct%coeffInd));
    !deallocate(residInd)
    
    !k = mexPrintf('Before O2'//achar(13))
    call findTempsForFit(O2Struct, TexStruct, coeff, dTCoeffs, T0Coeffs, Tex, dT0, T0);
    residInd = residIndEnd + (/(i, i = 1, size(O2Struct%data))/); residIndEnd = residInd(size(residInd))
    residual(residInd) = computeO2Residual(O2Struct, Tex, dT0, T0, coeff(O2Struct%coeffInd));
    !deallocate(residInd)

    !k = mexPrintf('Before rho'//achar(13))
    call findTempsForFit(rhoStruct, TexStruct, coeff, dTCoeffs, T0Coeffs, Tex, dT0, T0);
    residInd = residIndEnd + (/(i, i = 1, size(rhoStruct%data))/); residIndEnd = residInd(size(residInd))
    OlbDens = clamp(dble(10), evalMajorSpecies(rhoStruct, coeff(OStruct%coeffInd), OStruct%numBiases), &
                    dble(1E20));
    N2lbDens = clamp(dble(10), evalMajorSpecies(rhoStruct, coeff(N2Struct%coeffInd), N2Struct%numBiases), &
                     dble(1E20));
    HelbDens = clamp(dble(10), evalMajorSpecies(rhoStruct, coeff(HeStruct%coeffInd), HeStruct%numBiases), &
                     dble(1E20));
    ArlbDens = clamp(dble(10), evalMinorSpecies(rhoStruct, coeff(ArStruct%coeffInd), ArStruct%numBiases), &
                     dble(1E20));
    O2lbDens = clamp(dble(10.0), exp(coeff(O2Struct%coeffInd)), dble(1E20))
    modelRho = clamp(dble(1E-20), computeRho(T0, dT0, Tex, rhoStruct%Z, OlbDens, N2lbDens, HelbDens, &
                                             ArlbDens, O2lbDens), &
                     dble(0.1));
    residual(residInd) = (log(rhoStruct%data)/log(modelRho)) - 1;

    residual = weights * residual;

end function

function sumRowWise(vec, arr)
    implicit none
    real(kind = 8), intent(in) :: vec(:), arr(:,:)
    real(kind = 8), allocatable :: sumRowWise(:)
    integer :: i, numRows
    
    numRows = size(arr,1)
    allocate(sumRowWise(numRows))
    do i = 1, numRows
        sumRowWise(i) = dot_product(vec, arr(i,:))
    end do

end function

function clamp_scalar(minVal, vec, maxVal)
    real(kind = 8), intent(in) :: vec(:), minVal, maxVal
    real(kind = 8) :: clamp_scalar(size(vec))
    clamp_scalar = min(maxVal, max(minVal, vec))
end function

function clamp_vector(minVal, vec, maxVal)
    real(kind = 8), intent(in) :: vec(:), minVal(:), maxVal
    real(kind = 8) :: clamp_vector(size(vec))
    clamp_vector = min(maxVal, max(minVal, vec))
end function

end module fitModule

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use fitModule
    use lmSolverModule
    use omp_lib
    implicit none

    !     mexFunction arguments:
    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs

    !     Function declarations:
    mwPointer mxGetPr
    mwPointer mxCreateDoubleScalar
    mwPointer mxGetField
    mwPointer mxCreateDoubleMatrix
    mwPointer mxCreateString
    integer mxIsNumeric, mxIsStruct, mxIsDouble
    mwPointer mxGetM, mxGetN

    !     Pointers to input/output mxArrays:
    mwPointer x_ptr, y_ptr, field, rawArray,JTJ_ptr

    !     Array information:
    mwPointer mrows, ncols
    mwSize N, M

    real(kind = 8) :: tolX, tolFun, tolOpt, lambda0, firstOrderOpt
    integer :: maxFuncEvals, maxIter, exitFlag

    character(len = 200) :: tempChar

    integer(kind = 4) :: mexPrintf, k, mexEvalString, mexCallMATLAB, i,j

    real(kind = 8), allocatable :: y_output(:), solution(:), funVec(:), Jacobian(:,:), JTJ(:,:),JTJ_diag(:)

    !-----------------------------------------------------------------------
    !     Check for proper number of argumentS% 
    if(nrhs /= 11) then
        call mexErrMsgTxt ('Levenberg-Marquardt: Nine inputs required!')
    elseif(nlhs .gt. 2) then
        call mexErrMsgTxt ('Levenberg-Marquardt: At most two outputs allowed!')
    endif

    !     Validate inputs
    !     Check that the input is a struct.
    if(mxIsStruct(prhs(1)) .eq. 0 .or. mxIsStruct(prhs(2)) .eq. 0 .or. &
       mxIsStruct(prhs(3)) .eq. 0 .or. mxIsStruct(prhs(4)) .eq. 0 .or. &
       mxIsStruct(prhs(5)) .eq. 0 .or. mxIsStruct(prhs(6)) .eq. 0 .or. &
       mxIsStruct(prhs(7)) .eq. 0) then
        call mexErrMsgTxt ('Levenberg-Marquardt: Seven first inputs must be structs!')
    endif
    if (mxIsDouble(prhs(8)) .eq. 0 .or. mxIsDouble(prhs(9)) .eq. 0) then
        
    end if

    ! Read structs
    !k = mexPrintf('Before Tex'//achar(13))
    TexStruct = structToDerived_TexAndMajor(prhs(1), 'Tex')
    !k = mexPrintf('Before O'//achar(13))
    OStruct = structToDerived_TexAndMajor(prhs(2), 'O')
    !k = mexPrintf('Before N2'//achar(13))
    N2Struct = structToDerived_TexAndMajor(prhs(3), 'N2')
    !k = mexPrintf('Before He'//achar(13))
    HeStruct = structToDerived_TexAndMajor(prhs(4), 'He')
    !k = mexPrintf('Before Ar'//achar(13))
    ArStruct = structToDerived_TexAndMajor(prhs(5), 'Ar')
    !k = mexPrintf('Before O2'//achar(13))
    O2Struct = structToDerived_TexAndMajor(prhs(6), 'O2')
    !k = mexPrintf('Before Rho'//achar(13))
    rhoStruct = structToDerived_TexAndMajor(prhs(7), 'rho')
    
    !k = mexPrintf('Before dTCoeffs'//achar(13))
    N = max(mxGetM(prhs(8)), mxGetN(prhs(8)))
    allocate(dTCoeffs(N))
    call mxCopyPtrToReal8(mxGetPr(prhs(8)), dTCoeffs, N)
    !write(tempChar,*) dTCoeffs(N)
    !k = mexPrintf('last dT coeff: '//tempChar//achar(13))

    !k = mexPrintf('Before T0Coeffs'//achar(13))
    N = max(mxGetM(prhs(9)), mxGetN(prhs(9)))
    allocate(T0Coeffs(N))
    call mxCopyPtrToReal8(mxGetPr(prhs(9)), T0Coeffs, N)
    
    
    ! Read weights
    !k = mexPrintf('Before Weights'//achar(13))
    N = mxGetM(prhs(10))
    allocate(weights(N))
    call mxCopyPtrToReal8(mxGetPr(prhs(10)), weights, N)
    
    ! Read initGuess
    !k = mexPrintf('Before InitGuess'//achar(13))
    N = max(mxGetM(prhs(11)), mxGetN(prhs(11)))
    allocate(initGuess(N))
    call mxCopyPtrToReal8(mxGetPr(prhs(11)), initGuess, N)

    allocate(JTJ_diag(N))
    
    !$omp parallel
    if(omp_get_thread_num() == 0) write(tempChar,*) omp_get_num_threads()
    !$omp end parallel
    k = mexPrintf('OpenMP threads: '//tempChar//achar(13))
    
    
    ! ----- CALL Levenberg-Marquardt solver ------------
    tolX = 1E-8
    tolFun = 1E-5
    tolOpt = 1E4
    lambda0 = 1E-1
    maxFuncEvals = 5000 * size(initGuess)
    maxIter = 0 !!!!!!!!!!!!!!!!!!!!!!
    
    call lmSolve(modelMinimizationFunction, initGuess, tolX, tolFun, tolOpt, lambda0, maxFuncEvals, maxIter, &
                 JacobianAtSolution = Jacobian, solution = solution, funVec = funVec, exitFlag = exitFlag,&
                 firstOrderOptAtSolution = firstOrderOpt, JTJ_diag_sol = JTJ_diag)




    

    !k = mexPrintf('Before output'//achar(13))
    y_output = solution ! !!!!!!!!   
    !     Create matrix for the return argument.
    plhs(1) = mxCreateDoubleMatrix(size(y_output),1,0)
    plhs(2) = mxCreateDoubleMatrix(size(JTJ_diag),1,0)
    !k = mexPrintf('Before mxGetPr'//achar(13)) 
    y_ptr = mxGetPr(plhs(1))
    JTJ_ptr = mxGetPr(plhs(2))

    !     Load the data into y_ptr, which is the output to MATLAB.
    !k = mexPrintf('Before copy'//achar(13))
    !write(tempChar,*) size(y_output), mxGetM(plhs(1))
    !k = mexPrintf('Output size: '//tempChar//achar(13)) 
    
    call mxCopyReal8ToPtr(y_output,y_ptr,size(y_output))
    call mxCopyReal8ToPtr(JTJ_diag,JTJ_ptr,size(JTJ_diag))
    
    !k = mexPrintf('Begin deallocation'//achar(13))
    deallocate(weights, initGuess)
    call deallocateStruct(TexStruct, 'Tex')
    call deallocateStruct(OStruct, 'O')
    call deallocateStruct(N2Struct, 'N2')
    call deallocateStruct(HeStruct, 'He')
    call deallocateStruct(ArStruct, 'Ar')
    call deallocateStruct(O2Struct, 'O2')
    call deallocateStruct(rhoStruct, 'rho')
    k = mexPrintf('Deallocation performed'//achar(13))
    return
end

