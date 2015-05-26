!> @file
!>  DTM-2012 subroutines
!######################################################
!

!> This module provides the "model" custom type
module dtm_model

! version :  2012

  integer, parameter :: Nb_lat         =6
  integer, parameter :: Nb_flux        =12
  integer, parameter :: Nb_kp          =15
  integer, parameter :: Nb_SLat        =2
  integer, parameter :: Nb_SASLat      =2
  integer, parameter :: Nb_NSLat       =3
  integer, parameter :: Nb_SANSLat     =3
  integer, parameter :: Nb_DiAn        =12
  integer, parameter :: Nb_SDiAn       =12
  integer, parameter :: Nb_TDi         =2
  integer, parameter :: Nb_AMg         =9
  integer, parameter :: Nb_Lon         =10
  integer, parameter :: Nb_dPhas       =4

!> @brief
!> This type is used internally within the DTM2012 subroutine
  TYPE MODEL

                                                           !  various termes
   REAL(KIND(1.e0))                       :: T_one        ! SUM(T_x*T_dx) except T_dPhas
   REAL(KIND(1.e0)), dimension(Nb_lat)    :: T_Lat           !  in latitude
   REAL(KIND(1.e0)), dimension(Nb_flux)   :: T_flux         !  in flux
   REAL(KIND(1.e0)), dimension(Nb_kp)     :: T_kp         ! in kp
   REAL(KIND(1.e0)), dimension(Nb_SLat)   :: T_SLat         !  annual & symetric in latitude
   REAL(KIND(1.e0)), dimension(Nb_SASLat) :: T_SASLat         !  semi-annual & symetric in latitude
   REAL(KIND(1.e0)), dimension(Nb_NSLat)  :: T_NSLat         !  annual & non-symetric in latitude
   REAL(KIND(1.e0)), dimension(Nb_SANSLat):: T_SANSLat         !  semi-annual  non-symetric in latitude
   REAL(KIND(1.e0)), dimension(Nb_DiAn)   :: T_DiAn         !  diurnal (& annual coupled)
   REAL(KIND(1.e0)), dimension(Nb_SDiAn)  :: T_SDiAn         !  semi-diurnal (& annual coupled)
   REAL(KIND(1.e0)), dimension(Nb_TDi)    :: T_TDi         !  ter-diurnal
   REAL(KIND(1.e0)), dimension(Nb_AMg)    :: T_AMg         !  activity (Magnetic)
   REAL(KIND(1.e0)), dimension(Nb_Lon)    :: T_Lon         !  in longitude
   REAL(KIND(1.e0)), dimension(Nb_dPhas)  :: T_dPhas         !  derivative of phases of the annual and semi-annual terms


  END TYPE MODEL
 
    
end module dtm_model

!########################################################

!> This module provides variables shared among all DTM2012 subroutines
module dtm_model_par

 use  dtm_model

   !.. Structure .. 
   ! New DTM format/structure

  TYPE (MODEL)  ::  tt,h,he,ox,az2,o2,az,t0,tp   
  TYPE (MODEL)  ::  dtt,dh,dhe,dox,daz2,do2,daz,dt0,dtp  
  TYPE (MODEL)  :: MODEL_NULL


end module dtm_model_par


!!!!>>>  subroutine  dtm2012

 !---------------------------------------------------------------------------
 !
 ! ROUTINE: dtm2012
 !
 !> @author Sean Bruinsma
 !>
 !> @brief  calculation of temperature and density with DTM2012
 !
 ! PROTOTYPE:
 !                     dtm2012 (
 !                              real, intent(in) :: day
 !                              real, intent(in) :: f(2)
 !                              real, intent(in) :: fbar(2)
 !                              real, intent(in) :: akp(4)
 !                              real, intent(in) :: alti
 !                              real, intent(in) :: hl
 !                              real, intent(in) :: alat
 !                              real, intent(in) :: xlon
 !                              real, intent(out) :: tz
 !                              real, intent(out) :: tinf
 !                              real, intent(out) :: tp120
 !                              real, intent(out) :: ro
 !                              real, intent(out) :: d(6)
 !                              real, intent(out) :: wmm
 !                             )
 !
 !
 ! INPUT ARGUMENTS:
 !>                    @param[in] day day-of-year [1.-366.]
 !>                    @param[in] f(2) f(1)=instantaneous flux at (t - 24hr) / f(2)=0
 !>                    @param[in] fbar(2) fbar(1)=average flux of last 81 days at t  / fbar(2)=0.
 !>                    @param[in] akp(1): kp delayed by 3 hours
 !>                    @param[in] akp(2): must be 0.
 !>                    @param[in] akp(3): mean kp of last 24 hours 
 !>                    @param[in] akp(4): must be 0.
 !>                    @param[in] alti altitude (in km) greater than 120 km
 !>                    @param[in] hl local time (in radian: 0-24hr = 0-2pi)
 !>                    @param[in] alat latitude (in radian)
 !>                    @param[in] xlon longitude (in radian)
 !
 !
 ! OUTPUT ARGUMENTS:
 !>                    @param[out] ro           density in g/cm^3 at the given position
 !>                    @param[out] tinf         exosphere temperature
 !>                    @param[out] tz           temperature at the given height
 !>                    @param[out] xmm          mean molecular mass in gram
 !>                    @param[out] tp120        temperature gradient at 120 km (NOT accurate) 
 !>                    @param[out] d(6)         partial density of atomic hydrogen in g/cm^3 (1)
 !>                                             partial density of helium in g/cm^3 (2)
 !>                                             partial density of atomic oxygen in g/cm^3 (3)
 !>                                             partial density of molecular nitrogen in g/cm^3 (4)
 !>                                             partial density of molecular oxygen in g/cm^3 (5)
 !>                                             (6) not used
 !
 !> @date 06/2012
 !
 !---------------------------------------------------------------------------

subroutine dtm2012(day,f,fbar,akp,alti,hl,alat,xlon,tz,tinf,tp120,ro,d,wmm)

  !***********************************************************************
  !*aut SB         Mod PS     
  !*ver 06/10/2009 :: 06/2012
  !
  !*rol calculation of temperature and density with DTM2009/2012/2013  !*par ** INPUT **
  !*par ** INPUT **
  !     day      = day-of-year [1.-366.]
  !     f        = f(1)=instantaneous flux at (t - 24hr) / f(2)=0.
  !     fbar     = fbar(1)=average flux at t / fbar(2)=0.
  !     akp      = akp(1)= kp delayed by 3 hours, akp(3)=mean of last 24 hours / akp(2) & akp(4)=0.
  !     alti     = altitude (in km) greater than 120 km
  !     hl=      = local time (in radian: 0-24hr = 0-2pi)
  !     alat     = latitude (in radian)
  !     xlon     = longitude (in radian)

  !*par ** OUTPUT **
  !     tz        = temperature at altitude -> alti
  !     tinf      = exospheric temperature
  !     d(1)      = partial density of atomic hydrogen in g/cm3
  !     d(2)      = partial density of helium in g/cm3
  !     d(3)      = partial density of atomic oxygen in g/cm3
  !     d(4)      = partial density of molecular nitrogen in g/cm3
  !     d(5)      = partial density of molecular oxygen in g/cm3
  !     d(6)      = NOT USED
  !     ro       = density (total) in g/cm3
  !     wmm      = mean molecular mass
  !***********************************************************************
  !.. the uses .. 

    use  dtm_model_par

  ! 
  !.. Implicit Declarations .. 

  implicit none
  ! 
  !.. Parameters .. 
  ! 
  !.. Formal Arguments .. 
  real, intent(in) :: alat,alti,day,hl,xlon
  
  real, intent(out) :: ro,tinf,tz,tp120,wmm 

  real, dimension(2), intent(in) :: f,fbar
  real, dimension(4), intent(in) :: akp
  real, intent(out) :: d(6) 
  
  ! 
  !.. Local Scalars .. 
  integer :: i,i1,ialt
  integer :: ityp = 0
  integer :: i2,n
  real :: akp2,asc,attenuat,c,c2,c4,cecl4,cl,cmg,cp,cts,delkm,dt120,dtinf, &
          dtp120,gamma,gdelaz,gdelh,gdelo,gdelt0,gdeltp,s2,sasc,tinftz,ts,zeta
  real :: cpmg = .19081,re = 6356.77,rgas = 831.4

  real :: akp4,casc,cecl2,clmlmg,cmg2,cmg4,cpcl,cpsl,ctl,ctsg,difkm,dzeta, &
          dzeta2,expsz,gdelaz2,gdelhe,gdelo2,gdelt,glb,s,secl2,sigma,sigzeta, &
          sl,sp,stl,sts,stsg,t,t120sd,t120tt,t120tz,tp120sd,tsg,unsure,upapg

  real :: xlmg = -1.2392
  real :: t120 = 0.0,xlog = 0.0
  real :: cose = .9175,gsurf = 980.665,sine = .3978,spmg = .98163,zero = 0.
  real :: zlb
  real :: zlb0 = 120.
  ! 
  !.. Local Arrays .. 
  integer, dimension(6) :: ma
  real, dimension(6) :: alefa
  real, dimension(6) :: cc = 0.0,dbase = 0.0,fz = 0.0
  real, dimension(6) :: vma


!!*********************************************
  ! 
  !.. External Calls .. 
  external gldtm_XX

  !radg removed 18/09/2012. It is not used anywhere, and ifort complains about it
  !external tramkm

  ! 
  !.. Intrinsic Functions .. 
  intrinsic LOG, REAL, atan, cos, exp, int, mod, sin, sqrt
  ! 
  !.. Common Blocks .. 
  common /cons/ pi,deupi,cdr,sard
  ! 
  !... Variables in Common Block /cons/ ... 
  real :: cdr,deupi,pi,sard
  
  common /datmo/ npara,itype,ilin

  ! 
  !... Variables in Common Block /datmo/ ... 
  integer :: npara
  integer :: ilin,itype
  

  common /eclipt/ cecl,secl,c2ecl,s2ecl,c3ecl,s3ecl,p10ecl,p20ecl,p30ecl, &
                  p40ecl,p50ecl,p11ecl,p21ecl,p31ecl,p41ecl,p51ecl,p22ecl, &
                  p32ecl,p42ecl,p52ecl,p33ecl
  ! 
  !... Variables in Common Block /eclipt/ ... 
  real :: c2ecl,c3ecl,cecl,p10ecl,p11ecl,p20ecl,p21ecl,p22ecl,p30ecl,p31ecl, &
          p32ecl,p33ecl,p40ecl,p41ecl,p42ecl,p50ecl,p51ecl,p52ecl,s2ecl,s3ecl, &
          secl

  common /hlocal/ hl0,ch,sh,c2h,s2h,c3h,s3h
  ! 
  !... Variables in Common Block /hlocal/ ... 
  real :: c2h,c3h,ch,hl0,s2h,s3h,sh



!!*********************************************
  ! New DTM format/structure

  !TYPE (MODEL)  ::  tt,h,he,ox,az2,o2,az,t0,tp   
  !TYPE (MODEL)  ::  dtt,dh,dhe,dox,daz2,do2,daz,dt0,dtp  

  !common/pardtm/tt,h,he,ox,az2,o2,az,t0,tp

!!*********************************************

  common /plgdtm/ p10,p20,p30,p40,p50,p60,p11,p21,p31,p41,p51,p22,p32,p42,p52, &
                  p62,p33,p10mg,p20mg,p40mg
  ! 
  !... Variables in Common Block /plgdtm/ ... 
  real :: p10,p10mg,p11,p20,p20mg,p21,p22,p30,p31,p32,p33,p40,p40mg,p41,p42, &
          p50,p51,p52,p60,p62
  ! 
  !.. Data Declarations .. 
  data alefa/-0.40,-0.38,0.,0.,0.,0./
  data ma/1,4,16,28,32,14/
  data vma/ &
       1.6606e-24,6.6423e-24,26.569e-24,46.4958e-24,53.1381e-24,23.2479e-24/
  ! 
  ! ... Executable Statements ...
  ! 
  zlb = zlb0 
  !
  ialt = int(alti)
  ro = 0.

  dtinf = 0.
  dt120 = 0.
  dtp120 = 0.
  fz(1) = 0.
  fz(2) = 0.
  fz(3) = 0.
  fz(4) = 0.
  fz(5) = 0.
  fz(6) = 0.
  !
    dtt%T_one  = 0.
    dh%T_one   = 0.
    dhe%T_one  = 0.
    dox%T_one  = 0.
    daz2%T_one = 0.
    do2%T_one  = 0.
    daz%T_one  = 0.
    dt0%T_one  = 0.
    dtp%T_one  = 0.
  !
  !
  !   calcul des polynomes de legendre
  c = sin(alat)
  c2 = c * c
  c4 = c2 * c2
  s = cos(alat)
  s2 = s * s
  p10 = c
  p20 = 1.5*c2 - 0.5
  p30 = c * (2.5*c2-1.5)
  p40 = 4.375*c4 - 3.75*c2 + 0.375
  p50 = c * (7.875*c4-8.75*c2+1.875)
  p60 = (5.5*c*p50-2.5*p40) / 3.
  p11 = s
  p21 = 3. * c * s
  p31 = s * (7.5*c2-1.5)
  p41 = c * s * (17.5*c2-7.5)
  p51 = s * (39.375*c4-26.25*c2+1.875)
  p22 = 3. * s2
  p32 = 15. * c * s2
  p42 = s2 * (52.5*c2-7.5)
  p52 = 3.*c*p42 - 2.*p32
  p62 = 2.75*c*p52 - 1.75*p42
  p33 = 15. * s * s2
    !   calcul des polynomes de legendre / pole magnetique (79n,71w)
    clmlmg = cos(xlon-xlmg)
    sp = s*cpmg*clmlmg + c*spmg
    !
    cmg = sp ! pole magnetique
    cmg2 = cmg * cmg
    cmg4 = cmg2 * cmg2
    p10mg = cmg
    p20mg = 1.5*cmg2 - 0.5
    p40mg = 4.375*cmg4 - 3.75*cmg2 + 0.375
  !
  !   heure locale 
  hl0 = hl
  ch = cos(hl0)
  sh = sin(hl0)
  c2h = ch*ch - sh*sh
  s2h = 2. * ch * sh
  c3h = c2h*ch - s2h*sh
  s3h = s2h*ch + c2h*sh
  !
  !   calcul of fonction g(l) / tinf, t120, tp120
  !! new version declaration of tt,dtt en struct dtm
  call gldtm_XX(f,fbar,akp,day,tt,dtt,gdelt,1.,xlon)
  dtt%T_One = 1. + gdelt
  tinf = tt%T_One * dtt%T_One

  call gldtm_XX(f,fbar,akp,day,t0,dt0,gdelt0,1.,xlon)
  dt0%T_One = 1. + gdelt0
  t120 = t0%T_One * dt0%T_One

  call gldtm_XX(f,fbar,akp,day,tp,dtp,gdeltp,1.,xlon)
  dtp%T_One = 1. + gdeltp
  tp120 = tp%T_One * dtp%T_One  
  !-----------------------------------------------------------------------------
  !
  !   calcul of concentrations n(z): H, HE, O, N2, O2, N
  sigma = tp120 / (tinf-t120)
  dzeta = (re+zlb) / (re+alti)
  zeta = (alti-zlb) * dzeta
  dzeta2 = dzeta * dzeta
  sigzeta = sigma * zeta
  expsz = exp(-sigzeta)
  tz = tinf - (tinf-t120)*expsz

      call gldtm_XX(f,fbar,akp,day,h,dh,gdelh,0.,xlon)
       dh%T_One = exp(gdelh)
     dbase(1) = h%T_One * dh%T_One

      call gldtm_XX(f,fbar,akp,day,he,dhe,gdelhe,0.,xlon)
       dhe%T_One = exp(gdelhe)
      dbase(2) =  he%T_One * dhe%T_One

      call gldtm_XX(f,fbar,akp,day,ox,dox,gdelo,1.,xlon)
       dox%T_One = exp(gdelo)
      dbase(3) = ox%T_One * dox%T_One
 
      call gldtm_XX(f,fbar,akp,day,az2,daz2,gdelaz2,1.,xlon)
       daz2%T_One = exp(gdelaz2)
      dbase(4) = az2%T_One * daz2%T_One

      call gldtm_XX(f,fbar,akp,day,o2,do2,gdelo2,1.,xlon)
       do2%T_One = exp(gdelo2)
      dbase(5) =  o2%T_One * do2%T_One

      call gldtm_XX(f,fbar,akp,day,az,daz,gdelaz,1.,xlon)
       daz%T_One = exp(gdelaz)
      dbase(6) = az%T_One * daz%T_One
    !
    glb = gsurf / (1.+zlb/re)**2
    glb = glb / (sigma*rgas*tinf)
    t120tz = t120 / tz
    tinftz = tinf / tz
    do i = 1,6

      gamma = ma(i) * glb
      upapg = 1. + alefa(i) + gamma
      fz(i) = t120tz**upapg * exp(-sigzeta*gamma)
      !   concentrations en H, HE, O, N2, O2, N
      cc(i) = dbase(i) * fz(i)
      !   densites en H, HE, O, N2, O2, N
      d(i) = cc(i) * vma(i)
      !
      !   densite totale
      ro = ro + d(i)
    end do

      !  average of atomic mass                              
      wmm=ro/(vma(1)*(cc(1)+cc(2)+cc(3)+cc(4)+cc(5)+cc(6)))

end subroutine dtm2012


 !---------------------------------------------------------------------------
 !
 ! ROUTINE: gldtm_XX
 !
 !> @author Sean Bruinsma
 !>
 !> @brief  calculation of function g(l)  for dtm2009 & dtm2012
 !
 ! PROTOTYPE:
 !                     gldtm_XX (
 !                              real, intent(in) :: day
 !                              real, intent(in) :: ff0
 !                              real, intent(in) :: lon
 !                              real, intent(out) :: gdel
 !                              real, dimension(2), intent(in) :: f,fbar
 !                              real, dimension(4), intent(in) :: akp
 !                             )
 !
 ! INPUT ARGUMENTS:
 !>                    @param[in] day
 !>                    @param[in] ff0 (1 for oxygen, nitrogen, helium, temperature, 0 for hydrogen)
 !>                    @param[in] lon (currently not used)
 !>                    @param[in] f f(1)=instantaneous flux at (day - 1) / f(2)=0
 !>                    @param[in] fbar fbar(1)=average flux a t and derivative / fbar(2)=0
 !>                    @param[in] akp akp(1)= kp delayed by 3 hours, akp(3)=mean of last 24 hours / akp(2) & akp(4)=0
 !
 ! OUTPUT ARGUMENTS:
 !>                    @param[out] gdel result for g(l)
 !
 !> @date 06/03/2012
 !
 !---------------------------------------------------------------------------

!***********************************************************************
subroutine gldtm_XX(f,fbar,akp,day,DTM_12,ddtm_12,gdel,ff0,xlon)
  !***********************************************************************
  ! SB 06/03/2009 MOD PS 04/2012
  ! Corrections & new variables
  !*rol calculation of function g(l)  for dtm2009 & dtm2012
  !
  !     DTM =  table of coefficients for g(l) or for each component
  !     ddtm = table of derivatives dg(l)/da
  !     ff0=1 : for oxygen, nitrogen, helium, temperature
  !     ff0=0 : hydrogen
  !
  !     gdel=result for g(l)
  !     xlon=longitude (not used) 
  !***********************************************************************
  !
   !.. the uses .. 
   use  dtm_model

   !.. Implicit Declarations .. 

     implicit none

  ! 
  !.. Parameters .. >> in dtm_model
  ! 
  !.. Formal Arguments .. 

  integer :: kle_eq=0
  real, intent(in) :: day,ff0,xlon
  real, intent(out) :: gdel
  real, dimension(2), intent(in) :: f,fbar
  real, dimension(4), intent(in) :: akp

   ! 
   ! 
  ! Here de DTM 
  TYPE (MODEL)                  ::  DTM_12         !! the Old a()
  TYPE (MODEL)                  ::  ddtm_12        !! the Old da()
  ! 
  ! 
  !.. Local Scalars .. 
  integer :: i
  integer :: ikp = 0
  integer :: ikpm
  real :: a74,a77,a78,a88,f1f,fp
  real :: c2fi = 0.0
  real :: a89,a90,a91,aada,clfl,cos2te,coste,dakp,dakpm,dkp,dkpm,f0,fp1, &
          rsin2te,slfl
  real :: rot = .017214206,rot2 = .034428412,roth = .261799387, &
          rots = 7.27220e-05
  real :: rsinte
  ! 
  !.. Local Arrays .. 
  real, dimension(2) :: fbm150 = 0.0,fmfb = 0.0
  ! 
  !.. Intrinsic Functions .. 
  intrinsic abs, cos, sin
  ! 
  !.. Common Blocks .. 
  common /datmo/ npara,itype,ilin
  ! 
  !... Variables in Common Block /datmo/ ... 
  integer :: npara

  integer :: ilin,itype
  
  common /hlocal/ hl,ch,sh,c2h,s2h,c3h,s3h
  !     For hl: Not Read, Not Written
  !     For ch to s3h: Read, Not Written
  ! 
  !... Variables in Common Block /hlocal/ ... 
  real :: c2h,c3h,ch,hl,s2h,s3h,sh

  common /plgdtm/ p10,p20,p30,p40,p50,p60,p11,p21,p31,p41,p51,p22,p32,p42,p52, &
                  p62,p33,p10mg,p20mg,p40mg
  ! 
  !... Variables in Common Block /plgdtm/ ... 
  real :: p10,p10mg,p11,p20,p20mg,p21,p22,p30,p31,p32,p33,p40,p40mg,p41,p42, &
          p50,p51,p52,p60,p62
  ! 
  ! ... Executable Statements ...
  ! 
  !   termes in latitude

   ddtm_12%T_lat(1) = p20
   ddtm_12%T_lat(2) = p40
   ddtm_12%T_lat(3) = p10
   ddtm_12%T_lat(4) = p30
   ddtm_12%T_lat(5) = p50
   ddtm_12%T_lat(6) = p60

  !   termes of flux

  fmfb(1) = f(1) - fbar(1)
  fmfb(2) = f(2) - fbar(2)
  fbm150(1) = fbar(1) - 150.
  fbm150(2) = fbar(2)

  ddtm_12%T_flux(1) = fmfb(1)
  ddtm_12%T_flux(3) = fbm150(1)


    ddtm_12%T_flux(1) = ddtm_12%T_flux(1) + DTM_12%T_flux(5)*fmfb(2)
    ddtm_12%T_flux(3) = ddtm_12%T_flux(3) + DTM_12%T_flux(6)*fbm150(2)
    
    ddtm_12%T_flux(5) = fmfb(2) * (DTM_12%T_flux(1)+2.*DTM_12%T_flux(2)*ddtm_12%T_flux(1) &
                      + DTM_12%T_flux(7)*p10+DTM_12%T_flux(8)*p20+DTM_12%T_flux(9)*p30)

    ddtm_12%T_flux(6) = fbm150(2) * (DTM_12%T_flux(3)+2.*DTM_12%T_flux(4)*ddtm_12%T_flux(3) &
                      + DTM_12%T_flux(10)*p10+DTM_12%T_flux(11)*p20+DTM_12%T_flux(12)*p30)

  ddtm_12%T_flux(2)  = ddtm_12%T_flux(1) * ddtm_12%T_flux(1)
  ddtm_12%T_flux(4)  = ddtm_12%T_flux(3) * ddtm_12%T_flux(3)
  ddtm_12%T_flux(7)  = ddtm_12%T_flux(1) * p10
  ddtm_12%T_flux(8)  = ddtm_12%T_flux(1) * p20
  ddtm_12%T_flux(9)  = ddtm_12%T_flux(1) * p30
  ddtm_12%T_flux(10) = ddtm_12%T_flux(3) * p20
  ddtm_12%T_flux(11) = ddtm_12%T_flux(3) * p30
  ddtm_12%T_flux(12) = ddtm_12%T_flux(3) * p40
  !   termes in kp
    c2fi = 1. - p10mg*p10mg

    dkp = akp(1) + (DTM_12%T_kp(5)+c2fi*DTM_12%T_kp(6))*akp(2)
    dakp = DTM_12%T_kp(1) + DTM_12%T_kp(2)*p20mg + DTM_12%T_kp(11)*p40mg + &
           2.*dkp*(DTM_12%T_kp(3)+DTM_12%T_kp(4)*p20mg+DTM_12%T_kp(14)*2.*dkp*dkp)

    ddtm_12%T_kp(5) = dakp * akp(2)
    ddtm_12%T_kp(6) =  ddtm_12%T_kp(5) * c2fi

    dkpm = akp(3) + DTM_12%T_kp(10)*akp(4)

    dakpm = DTM_12%T_kp(7) + DTM_12%T_kp(8)*p20mg + DTM_12%T_kp(12)*p40mg + &
            2.*dkpm*(DTM_12%T_kp(9)+DTM_12%T_kp(13)*p20mg + & 
          DTM_12%T_kp(15)*2.*dkpm*dkpm)

    ddtm_12%T_kp(10) = dakpm * akp(4)

  ddtm_12%T_kp(1) = dkp
  ddtm_12%T_kp(2) = p20mg * dkp
  ddtm_12%T_kp(11) = p40mg * dkp
  ddtm_12%T_kp(3) = dkp * dkp
  ddtm_12%T_kp(4) = p20mg * ddtm_12%T_kp(3)
  ddtm_12%T_kp(14) = ddtm_12%T_kp(3) * ddtm_12%T_kp(3)
  ddtm_12%T_kp(7) = dkpm
  ddtm_12%T_kp(8) = p20mg * dkpm
  ddtm_12%T_kp(12) = p40mg * dkpm
  ddtm_12%T_kp(9) = dkpm * dkpm
  ddtm_12%T_kp(13) = p20mg * ddtm_12%T_kp(9)
  ddtm_12%T_kp(15) = ddtm_12%T_kp(9) * ddtm_12%T_kp(9)

  !   function: g(l) non periodic

  f0 = DTM_12%T_flux(1)*ddtm_12%T_flux(1) & 
     + DTM_12%T_flux(2)*ddtm_12%T_flux(2) & 
     + DTM_12%T_flux(3)*ddtm_12%T_flux(3) &
     + DTM_12%T_flux(4)*ddtm_12%T_flux(4) &
     + DTM_12%T_flux(7)*ddtm_12%T_flux(7)  &
     + DTM_12%T_flux(8)*ddtm_12%T_flux(8)  &
     + DTM_12%T_flux(9)*ddtm_12%T_flux(9)  &
     + DTM_12%T_flux(10)*ddtm_12%T_flux(10)  &
     + DTM_12%T_flux(11)*ddtm_12%T_flux(11)  &
     + DTM_12%T_flux(12)*ddtm_12%T_flux(12)

  f1f = 1. + f0*ff0

  !     write(6,*) ' ff0,f0,f1f ',ff0,f0,f1f

  f0 = f0 + DTM_12%T_lat(1)*ddtm_12%T_lat(1) &
          + DTM_12%T_lat(2)*ddtm_12%T_lat(2) &
        + DTM_12%T_lat(3)*ddtm_12%T_lat(3) &
        + DTM_12%T_lat(4)*ddtm_12%T_lat(4) &
        + DTM_12%T_kp(1)*ddtm_12%T_kp(1) &
        + DTM_12%T_kp(2)*ddtm_12%T_kp(2) &
        + DTM_12%T_kp(3)*ddtm_12%T_kp(3) &
        + DTM_12%T_kp(4)*ddtm_12%T_kp(4) &
        + DTM_12%T_kp(11)*ddtm_12%T_kp(11) &
        + DTM_12%T_kp(7)*ddtm_12%T_kp(7) &
          + DTM_12%T_kp(8)*ddtm_12%T_kp(8) &
        + DTM_12%T_kp(9)*ddtm_12%T_kp(9) &
        + DTM_12%T_kp(12)*ddtm_12%T_kp(12) &
        + DTM_12%T_kp(13)*ddtm_12%T_kp(13) &
          + DTM_12%T_kp(14)*ddtm_12%T_kp(14) & 
        + DTM_12%T_kp(15) *ddtm_12%T_kp(15)  &
        + DTM_12%T_lat(5)*ddtm_12%T_lat(5) &
        + DTM_12%T_lat(6)*ddtm_12%T_lat(6)

  !   terms annual & symetric in latitude 
  ddtm_12%T_SLat(1) = cos(rot*(day-DTM_12%T_dPhas(1)))
  ddtm_12%T_SLat(2) = p20 * ddtm_12%T_SLat(1)
  !   terms semi-annual & symetric in latitude
  ddtm_12%T_SASLat(1) = cos(rot2*(day-DTM_12%T_dPhas(2)))
  ddtm_12%T_SASLat(2) = p20 * ddtm_12%T_SASLat(1)
  !   terms annual & non-symetric in latitude 
  coste = cos(rot*(day-DTM_12%T_dPhas(3)))
  ddtm_12%T_NSLat(1) = p10 * coste
  ddtm_12%T_NSLat(2) = p30 * coste
  ddtm_12%T_NSLat(3) = ddtm_12%T_flux(3) * ddtm_12%T_NSLat(1)
  !      ddtm_12%T_NSLat(3)=p50*coste
  !   terms semi-annual  non-symetric in latitude 
  cos2te = cos(rot2*(day-DTM_12%T_dPhas(4)))
  ddtm_12%T_SANSLat(1) = p10 * cos2te
  ddtm_12%T_SANSLat(2) = p30 * cos2te
  ddtm_12%T_SANSLat(3) = ddtm_12%T_flux(3) * ddtm_12%T_SANSLat(1)
  !      ddtm_12%T_SANSLat(3)=p50*cos2te
  !   terms diurnal (& annual coupled)
  ddtm_12%T_DiAn(1) = p11 * ch
  ddtm_12%T_DiAn(2) = p31 * ch
  ddtm_12%T_DiAn(3) = ddtm_12%T_flux(3) * ddtm_12%T_DiAn(1) 
  ddtm_12%T_DiAn(4) = ddtm_12%T_DiAn(1) * coste
  ddtm_12%T_DiAn(5) = p21 * ch * coste
  ddtm_12%T_DiAn(6) = p11 * sh
  ddtm_12%T_DiAn(7) = p31 * sh
  ddtm_12%T_DiAn(8) = ddtm_12%T_flux(3) * ddtm_12%T_DiAn(6)
  ddtm_12%T_DiAn(9) = ddtm_12%T_DiAn(6) * coste
  ddtm_12%T_DiAn(10) = p21 * sh * coste
  ddtm_12%T_DiAn(11)=p51*ch  
  ddtm_12%T_DiAn(12)=p51*sh  
  !   terms semi-diurnes (& annual coupled)
  ddtm_12%T_SDiAn(1) = p22 * c2h
  ddtm_12%T_SDiAn(2) = p42 * c2h
  ddtm_12%T_SDiAn(3) = p32 * c2h * coste
  ddtm_12%T_SDiAn(4) = p22 * s2h
  ddtm_12%T_SDiAn(5) = p42 * s2h
  ddtm_12%T_SDiAn(6) = p32 * s2h * coste
  ddtm_12%T_SDiAn(7) = p32 * c2h !coeff. rajoute pour tp120/t120 (slb)
  ddtm_12%T_SDiAn(8) = p32 * s2h
  ddtm_12%T_SDiAn(9) = ddtm_12%T_flux(3) * ddtm_12%T_SDiAn(1)
  ddtm_12%T_SDiAn(10) = ddtm_12%T_flux(3) * ddtm_12%T_SDiAn(4)
  ddtm_12%T_SDiAn(11) = p62 * c2h
  ddtm_12%T_SDiAn(12) = p62 * s2h
  !   terms ter-diurnes
  ddtm_12%T_TDi(1) = p33 * c3h
  ddtm_12%T_TDi(2) = p33 * s3h
  !   function periodic -> g(l) 
  fp = DTM_12%T_SLat(1)*ddtm_12%T_SLat(1) &
     + DTM_12%T_SLat(2)*ddtm_12%T_SLat(2) &
     + DTM_12%T_SASLat(1)*ddtm_12%T_SASLat(1) &
     + DTM_12%T_SASLat(2)*ddtm_12%T_SASLat(2) &
     + DTM_12%T_NSLat(1)*ddtm_12%T_NSLat(1)  &
     + DTM_12%T_NSLat(2)*ddtm_12%T_NSLat(2) & 
     + DTM_12%T_NSLat(3)*ddtm_12%T_NSLat(3)  &
     + DTM_12%T_SANSLat(1)*ddtm_12%T_SANSLat(1) &
     + DTM_12%T_DiAn(1)*ddtm_12%T_DiAn(1)   &
     + DTM_12%T_DiAn(2)*ddtm_12%T_DiAn(2) &
     + DTM_12%T_DiAn(3)*ddtm_12%T_DiAn(3)   &
     + DTM_12%T_DiAn(4)*ddtm_12%T_DiAn(4) &
     + DTM_12%T_DiAn(5)*ddtm_12%T_DiAn(5)   &
     + DTM_12%T_DiAn(6)*ddtm_12%T_DiAn(6) &
     + DTM_12%T_DiAn(7)*ddtm_12%T_DiAn(7)   &
     + DTM_12%T_DiAn(8)*ddtm_12%T_DiAn(8) &
     + DTM_12%T_DiAn(9)*ddtm_12%T_DiAn(9)   &
     + DTM_12%T_DiAn(10)*ddtm_12%T_DiAn(10) &
     + DTM_12%T_SDiAn(1)*ddtm_12%T_SDiAn(1)  &
     + DTM_12%T_SDiAn(3)*ddtm_12%T_SDiAn(3) &
     + DTM_12%T_SDiAn(4)*ddtm_12%T_SDiAn(4)  &
     + DTM_12%T_SDiAn(6) *ddtm_12%T_SDiAn(6) &
     + DTM_12%T_TDi(1)*ddtm_12%T_TDi(1)    &
     + DTM_12%T_TDi(2)*ddtm_12%T_TDi(2) &
     + DTM_12%T_SDiAn(2)*ddtm_12%T_SDiAn(2)  &
     + DTM_12%T_SDiAn(5)*ddtm_12%T_SDiAn(5) &
     + DTM_12%T_SANSLat(2)*ddtm_12%T_SANSLat(2) &
     + DTM_12%T_SANSLat(3)*ddtm_12%T_SANSLat(3) &
     + DTM_12%T_SDiAn(7)*ddtm_12%T_SDiAn(7)  &
     + DTM_12%T_SDiAn(8)*ddtm_12%T_SDiAn(8)  &
     + DTM_12%T_SDiAn(9)*ddtm_12%T_SDiAn(9)  &
     + DTM_12%T_SDiAn(10)*ddtm_12%T_SDiAn(10) &
     + DTM_12%T_SDiAn(11) *ddtm_12%T_SDiAn(11) &
     + DTM_12%T_SDiAn(12)*ddtm_12%T_SDiAn(12)  &
     + DTM_12%T_DiAn(11)*ddtm_12%T_DiAn(11)  &
     + DTM_12%T_DiAn(12)*ddtm_12%T_DiAn(12) 
  !
  !   terms magnetic activity 
  
  dDTM_12%T_AMg(1) = p10 * coste * dkp
  dDTM_12%T_AMg(2) = p30 * coste * dkp
  dDTM_12%T_AMg(3) = p50 * coste * dkp
  dDTM_12%T_AMg(4) = p11 * ch * dkp
  dDTM_12%T_AMg(5) = p31 * ch * dkp
  dDTM_12%T_AMg(6) = p51 * ch * dkp
  dDTM_12%T_AMg(7) = p11 * sh * dkp
  dDTM_12%T_AMg(8) = p31 * sh * dkp
  dDTM_12%T_AMg(9) = p51 * sh * dkp
  !
  !   function g(l) (additional periodic)
  fp = fp + DTM_12%T_AMg(1)*dDTM_12%T_AMg(1) &
          + DTM_12%T_AMg(2)*dDTM_12%T_AMg(2) &
        + DTM_12%T_AMg(3)*dDTM_12%T_AMg(3) &
        + DTM_12%T_AMg(4)*dDTM_12%T_AMg(4) &
          + DTM_12%T_AMg(5)*dDTM_12%T_AMg(5) &
        + DTM_12%T_AMg(6)*dDTM_12%T_AMg(6) &
        + DTM_12%T_AMg(7)*dDTM_12%T_AMg(7) &
        + DTM_12%T_AMg(8)*dDTM_12%T_AMg(8) &
          + DTM_12%T_AMg(9)*dDTM_12%T_AMg(9)
  !
    dakp = (DTM_12%T_AMg(1)*p10+DTM_12%T_AMg(2)*p30+DTM_12%T_AMg(3)*p50)*coste &
         + (DTM_12%T_AMg(4)*p11+DTM_12%T_AMg(5)*p31+DTM_12%T_AMg(6)*p51)*ch  &
         + (DTM_12%T_AMg(7)*p11+DTM_12%T_AMg(8)*p31+DTM_12%T_AMg(9)*p51)*sh

    dDTM_12%T_kp(5) = dDTM_12%T_kp(5) + dakp*akp(2)
    dDTM_12%T_kp(6) = dDTM_12%T_kp(5) + dakp*c2fi*akp(2)

    !   terms in longitude

    clfl = cos(xlon)
    dDTM_12%T_Lon(1) = p11 * clfl
    dDTM_12%T_Lon(2) = p21 * clfl
    dDTM_12%T_Lon(3) = p31 * clfl
    dDTM_12%T_Lon(4) = p41 * clfl
    dDTM_12%T_Lon(5) = p51 * clfl
    slfl = sin(xlon)
    dDTM_12%T_Lon(6) = p11 * slfl
    dDTM_12%T_Lon(7) = p21 * slfl
    dDTM_12%T_Lon(8) = p31 * slfl
    dDTM_12%T_Lon(9) = p41 * slfl
    dDTM_12%T_Lon(10) = p51 * slfl
    !
    !   function g(l) periodic (additional)
    Do i=1 , Nb_Lon
     fp = fp  + DTM_12%T_Lon(i)*dDTM_12%T_Lon(i)  
    EndDo
  !
  !   function g(l) sum (coupled with flux)

  gdel = f0 + fp*f1f
  
  !
  !
 ! --------------------------------------------------
end subroutine gldtm_XX


 !---------------------------------------------------------------------------
 !
 ! ROUTINE: P_ReadDTM12
 !
 !> @author Sean Bruinsma
 !>
 !> @brief  read DTM format/version 2012
 !
 ! PROTOTYPE:
 !                 P_ReadDTM12 (
 !                             )
 !
 !
 ! INPUT ARGUMENTS:
 !
 !
 ! OUTPUT ARGUMENTS:
 !
 !> @date 03/2012
 !
 !---------------------------------------------------------------------------
!***********************************************************************
subroutine P_ReadDTM12
  !***********************************************************************
  !
  !*rol Programme: read DTM format/version 2012 
  !
  !      pgf90  -o P_ReadDTM12.x P_ReadDTM12.f90
  !
  !*vers      : 03/2012
  !*aut      : PS
  !
!***********************************************************************
  !

  !
  !.. the uses .. 

   use  dtm_model_par
 
  !.. Implicit Declarations .. 
  implicit none


  !.. Parameters .. 
! _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.

  integer, parameter :: IVERIF = 0  !! IVERIF = 0/1 >> print model in u_Out
! _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.


  !.. VARIABLES .. 

  integer :: i,j,j_test,ik, npdtm,u_in
  integer :: Iok
  character titre*100,ct120*8,ctp120*8,Com*190

  
  ! New DTM format/structure

  !TYPE (MODEL)  ::  tt,h,he,ox,az2,o2,az,t0,tp   
  !TYPE (MODEL)  ::  dtt,dh,dhe,dox,daz2,do2,daz,dt0,dtp  

! --------- C o m m o n s -----------------------------------------


! --------------------------------------------------
  u_In=1
  open(u_In, file='dtm_2013_NF.dat',form = 'formatted',status='OLD')

! --------------------------------------------------
!***********  read 

      read(u_In,530) titre
!      write(*,530) titre
      read(u_In,*) npdtm

      backspace(u_In)

      read(u_In,546)  Com

! --------------------------------------------------

    Iok = 1    ! test read ok/nok
     
! ...... termes ONE
    
    i=1
    READ(u_In,640) j_test,  tt%T_one  ,dtt%T_one, &
                        h%T_one   ,dh%T_one, &
                  he%T_one  ,dhe%T_one, &
                  ox%T_one  ,dox%T_one, &
                  az2%T_one ,daz2%T_one, &
                  o2%T_one  ,do2%T_one, &
                  az%T_one  ,daz%T_one, &
                  t0%T_one  ,dt0%T_one, &
                  tp%T_one  ,dtp%T_one

!! ..:
 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_one " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_lat & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
   goto 15
 EndIf



! ...... termes in LAT
   Do i=1,Nb_lat

!!    READ(u_In,640) i,(DTM_12(j)%T_lat(i),ddtm_12(j)%T_lat(i),j=1,9)

    READ(u_In,640) j_test,  &
                       tt%T_lat(i)  ,dtt%T_lat(i), &
                        h%T_lat(i)   ,dh%T_lat(i), &
                  he%T_lat(i)  ,dhe%T_lat(i), &
                  ox%T_lat(i)  ,dox%T_lat(i), &
                  az2%T_lat(i) ,daz2%T_lat(i), &
                  o2%T_lat(i)  ,do2%T_lat(i), &
                  az%T_lat(i)  ,daz%T_lat(i), &
                  t0%T_lat(i)  ,dt0%T_lat(i), &
                  tp%T_lat(i)  ,dtp%T_lat(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_lat " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_lat & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo

! --------------------------------------------------
! ...... termes in Flux
   Do i=1,Nb_flux

    READ(u_In,640) j_test,  &
                       tt%T_flux(i)  ,dtt%T_flux(i), &
                        h%T_flux(i)   ,dh%T_flux(i), &
                  he%T_flux(i)  ,dhe%T_flux(i), &
                  ox%T_flux(i)  ,dox%T_flux(i), &
                  az2%T_flux(i) ,daz2%T_flux(i), &
                  o2%T_flux(i)  ,do2%T_flux(i), &
                  az%T_flux(i)  ,daz%T_flux(i), &
                  t0%T_flux(i)  ,dt0%T_flux(i), &
                  tp%T_flux(i)  ,dtp%T_flux(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_flux " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_flux & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo

! --------------------------------------------------
! ...... termes in kp
   Do i=1,Nb_kp

    READ(u_In,640) j_test,  &
                       tt%T_kp(i)  ,dtt%T_kp(i), &
                        h%T_kp(i)   ,dh%T_kp(i), &
                  he%T_kp(i)  ,dhe%T_kp(i), &
                  ox%T_kp(i)  ,dox%T_kp(i), &
                  az2%T_kp(i) ,daz2%T_kp(i), &
                  o2%T_kp(i)  ,do2%T_kp(i), &
                  az%T_kp(i)  ,daz%T_kp(i), &
                  t0%T_kp(i)  ,dt0%T_kp(i), &
                  tp%T_kp(i)  ,dtp%T_kp(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_kp " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_kp & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo

! --------------------------------------------------
! ...... termes in SLat
   Do i=1,Nb_SLat

    READ(u_In,640) j_test,  &
                       tt%T_SLat(i)  ,dtt%T_SLat(i), &
                        h%T_SLat(i)   ,dh%T_SLat(i), &
                  he%T_SLat(i)  ,dhe%T_SLat(i), &
                  ox%T_SLat(i)  ,dox%T_SLat(i), &
                  az2%T_SLat(i) ,daz2%T_SLat(i), &
                  o2%T_SLat(i)  ,do2%T_SLat(i), &
                  az%T_SLat(i)  ,daz%T_SLat(i), &
                  t0%T_SLat(i)  ,dt0%T_SLat(i), &
                  tp%T_SLat(i)  ,dtp%T_SLat(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_SLat " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_SLat & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo

! --------------------------------------------------

! ...... termes in SASLat
   Do i=1,Nb_SASLat

    READ(u_In,640) j_test,  &
                       tt%T_SASLat(i)  ,dtt%T_SASLat(i), &
                        h%T_SASLat(i)   ,dh%T_SASLat(i), &
                  he%T_SASLat(i)  ,dhe%T_SASLat(i), &
                  ox%T_SASLat(i)  ,dox%T_SASLat(i), &
                  az2%T_SASLat(i) ,daz2%T_SASLat(i), &
                  o2%T_SASLat(i)  ,do2%T_SASLat(i), &
                  az%T_SASLat(i)  ,daz%T_SASLat(i), &
                  t0%T_SASLat(i)  ,dt0%T_SASLat(i), &
                  tp%T_SASLat(i)  ,dtp%T_SASLat(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_SASLat " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_SASLat & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo

! --------------------------------------------------

! ...... termes in T_NSLat
   Do i=1,Nb_NSLat

    READ(u_In,640) j_test,  &
                       tt%T_NSLat(i)  ,dtt%T_NSLat(i), &
                        h%T_NSLat(i)   ,dh%T_NSLat(i), &
                  he%T_NSLat(i)  ,dhe%T_NSLat(i), &
                  ox%T_NSLat(i)  ,dox%T_NSLat(i), &
                  az2%T_NSLat(i) ,daz2%T_NSLat(i), &
                  o2%T_NSLat(i)  ,do2%T_NSLat(i), &
                  az%T_NSLat(i)  ,daz%T_NSLat(i), &
                  t0%T_NSLat(i)  ,dt0%T_NSLat(i), &
                  tp%T_NSLat(i)  ,dtp%T_NSLat(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_NSLat " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_NSLat & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo
! --------------------------------------------------

! ...... termes in T_SANSLat 
   Do i=1,Nb_SANSLat

    READ(u_In,640) j_test,  &
                       tt%T_SANSLat(i)  ,dtt%T_SANSLat(i), &
                        h%T_SANSLat(i)   ,dh%T_SANSLat(i), &
                  he%T_SANSLat(i)  ,dhe%T_SANSLat(i), &
                  ox%T_SANSLat(i)  ,dox%T_SANSLat(i), &
                  az2%T_SANSLat(i) ,daz2%T_SANSLat(i), &
                  o2%T_SANSLat(i)  ,do2%T_SANSLat(i), &
                  az%T_SANSLat(i)  ,daz%T_SANSLat(i), &
                  t0%T_SANSLat(i)  ,dt0%T_SANSLat(i), &
                  tp%T_SANSLat(i)  ,dtp%T_SANSLat(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_SANSLat " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_SANSLat & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo
! --------------------------------------------------
! ...... termes in DiAn
   Do i=1,Nb_DiAn

    READ(u_In,640) j_test,  &
                       tt%T_DiAn(i)  ,dtt%T_DiAn(i), &
                        h%T_DiAn(i)   ,dh%T_DiAn(i), &
                  he%T_DiAn(i)  ,dhe%T_DiAn(i), &
                  ox%T_DiAn(i)  ,dox%T_DiAn(i), &
                  az2%T_DiAn(i) ,daz2%T_DiAn(i), &
                  o2%T_DiAn(i)  ,do2%T_DiAn(i), &
                  az%T_DiAn(i)  ,daz%T_DiAn(i), &
                  t0%T_DiAn(i)  ,dt0%T_DiAn(i), &
                  tp%T_DiAn(i)  ,dtp%T_DiAn(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_DiAn " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_DiAn & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo
! --------------------------------------------------
! ...... termes  SDiAn
   Do i=1,Nb_SDiAn

    READ(u_In,640) j_test,  &
                       tt%T_SDiAn(i)  ,dtt%T_SDiAn(i), &
                        h%T_SDiAn(i)   ,dh%T_SDiAn(i), &
                  he%T_SDiAn(i)  ,dhe%T_SDiAn(i), &
                  ox%T_SDiAn(i)  ,dox%T_SDiAn(i), &
                  az2%T_SDiAn(i) ,daz2%T_SDiAn(i), &
                  o2%T_SDiAn(i)  ,do2%T_SDiAn(i), &
                  az%T_SDiAn(i)  ,daz%T_SDiAn(i), &
                  t0%T_SDiAn(i)  ,dt0%T_SDiAn(i), &
                  tp%T_SDiAn(i)  ,dtp%T_SDiAn(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_SDiAn " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_SDiAn & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo
! ------------- -------------------------------------
! ...... termes in  TDi
   Do i=1,Nb_TDi

    READ(u_In,640) j_test,  &
                       tt%T_TDi(i)  ,dtt%T_TDi(i), &
                        h%T_TDi(i)   ,dh%T_TDi(i), &
                  he%T_TDi(i)  ,dhe%T_TDi(i), &
                  ox%T_TDi(i)  ,dox%T_TDi(i), &
                  az2%T_TDi(i) ,daz2%T_TDi(i), &
                  o2%T_TDi(i)  ,do2%T_TDi(i), &
                  az%T_TDi(i)  ,daz%T_TDi(i), &
                  t0%T_TDi(i)  ,dt0%T_TDi(i), &
                  tp%T_TDi(i)  ,dtp%T_TDi(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_TDi " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_TDi & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo
! --------------------------------------------------
! ...... termes in  AMg
   Do i=1,Nb_AMg

    READ(u_In,640) j_test,  &
                       tt%T_AMg(i)  ,dtt%T_AMg(i), &
                        h%T_AMg(i)   ,dh%T_AMg(i), &
                  he%T_AMg(i)  ,dhe%T_AMg(i), &
                  ox%T_AMg(i)  ,dox%T_AMg(i), &
                  az2%T_AMg(i) ,daz2%T_AMg(i), &
                  o2%T_AMg(i)  ,do2%T_AMg(i), &
                  az%T_AMg(i)  ,daz%T_AMg(i), &
                  t0%T_AMg(i)  ,dt0%T_AMg(i), &
                  tp%T_AMg(i)  ,dtp%T_AMg(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_AMg " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_AMg & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo
! --------------------------------------------------
! ...... termes in  Lon
   Do i=1,Nb_Lon

    READ(u_In,640) j_test,  &
                       tt%T_Lon(i)  ,dtt%T_Lon(i), &
                        h%T_Lon(i)   ,dh%T_Lon(i), &
                  he%T_Lon(i)  ,dhe%T_Lon(i), &
                  ox%T_Lon(i)  ,dox%T_Lon(i), &
                  az2%T_Lon(i) ,daz2%T_Lon(i), &
                  o2%T_Lon(i)  ,do2%T_Lon(i), &
                  az%T_Lon(i)  ,daz%T_Lon(i), &
                  t0%T_Lon(i)  ,dt0%T_Lon(i), &
                  tp%T_Lon(i)  ,dtp%T_Lon(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_Lon " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_Lon & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo
!
! ...... termes in dPhas 
   Do i=1,Nb_dPhas

    READ(u_In,640) j_test,  &
                       tt%T_dPhas(i)  ,dtt%T_dPhas(i), &
                        h%T_dPhas(i)   ,dh%T_dPhas(i), &
                  he%T_dPhas(i)  ,dhe%T_dPhas(i), &
                  ox%T_dPhas(i)  ,dox%T_dPhas(i), &
                  az2%T_dPhas(i) ,daz2%T_dPhas(i), &
                  o2%T_dPhas(i)  ,do2%T_dPhas(i), &
                  az%T_dPhas(i)  ,daz%T_dPhas(i), &
                  t0%T_dPhas(i)  ,dt0%T_dPhas(i), &
                  tp%T_dPhas(i)  ,dtp%T_dPhas(i)

 If (i .NE.  j_test) THEN
  Write(*,*) "***WARNING: PB in field T_dPhas " 
  Write(*,*) "***WARNING: incompatibility beetween dimension T_dPhas & DTM file  " 
  Write(*,*) "*** i=",i,' j_test=', j_test
     Iok = 0 
  goto 15
 EndIf

    EndDo

! --------------------------------------------------

  530   format(a100)
  546   format(a190)

! --------------------------------------------------

close(u_In)

15   continue

      If (Iok.EQ.0) THEN
        stop
       write(*,*) " exit: LecDTM"
      ENDIF

! --------------------------------------------------
 640   format(i4,9(e13.6,e9.2))
 
! --------------------------------------------------
end subroutine P_ReadDTM12
!
!***********************************************************************
