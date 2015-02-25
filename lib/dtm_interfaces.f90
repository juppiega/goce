!> @file
!> Header definitions for DTM2013 wrappers


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
!> @brief           Explicit interfaces for dtm_wrapper
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
! Software context            : ATMOP                                        !
! Library Responsible         : Noelia Sanchez                               !
!>Subroutine Author           : @author Raul Dominguez
! Company                     : DEIMOS Space S.L.                            !
! Programming Language        : Fortran 90                                   !
! Associated File Name        : dtm_wrapper.f90                              !
! Development Compiler and OS : Windows 7 32 bit, cygwin                     !
! Compiler Version            : GNU gfortran 4.5.3                           !
! Compiling Options           : standard                                     !
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
!  Method                                                                    !
!  ======
!>  @details
!>   This module provides explicit interfaces for the dtm_wrapper subroutine
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
! Function history                                                           !
! ================                                                           !
!                                                                            !
!> @version 1.0
!> @date 15/07/2013
!  Initial coding
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
module dtm_wrapper_interfaces

    use t_dtm_date
    implicit none

    private
    public dtm_wrapper_interface
    public dtm_wrapper
    public load_config_interface
    public load_config

    interface dtm_wrapper_interface
     subroutine dtm_wrapper(in_date,alti,alat,xlon,tz,tinf,tp120,ro,ro_unc,d,wmm, &
                           f_out,fbar_out,akp_out,hl_out,dayofyear_out)
       use t_dtm_date
       type(dtm_date), intent(inout) :: in_date
       real, intent(in) :: alat
       real, intent(in) :: alti
       real, intent(in) :: xlon
       real, intent(out) :: ro
       real, intent(out) :: ro_unc
       real, intent(out) :: tinf
       real, intent(out) :: tz
       real, intent(out) :: tp120
       real, intent(out) :: wmm
       real, intent(out), dimension(6) :: d
       real, dimension(2), optional, intent(out) :: f_out
       real, dimension(2), optional, intent(out) :: fbar_out
       real, dimension(4), optional, intent(out) :: akp_out
       real, optional, intent(out) :: hl_out
       real, optional, intent(out) :: dayofyear_out
     end subroutine
    end interface

    interface load_config_interface
      subroutine load_config(unit,file)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: file
      end subroutine
    end interface

end module dtm_wrapper_interfaces


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
!> @brief           Explicit interfaces for main DTM2013 routines
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
! Software context            : ATMOP                                        !
! Library Responsible         : Noelia Sanchez                               !
!>Subroutine Author           : @author Raul Dominguez
! Company                     : DEIMOS Space S.L.                            !
! Programming Language        : Fortran 90                                   !
! Associated File Name        : dtm_wrapper.f90                              !
! Development Compiler and OS : Windows 7 32 bit, cygwin                     !
! Compiler Version            : GNU gfortran 4.5.3                           !
! Compiling Options           : standard                                     !
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
!  Method                                                                    !
!  ======
!>  @details
!>   This module provides explicit interfaces for the dtm2013,
!>   density_uncertainty, and P_ReadDTM12 subroutines
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
! Function history                                                           !
! ================                                                           !
!                                                                            !
!> @version 1.0
!> @date 15/11/2013
!  Upgraded from DTM2012
!                                                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
module dtm2013_interfaces

    implicit none

    private
    public dtm2013_interface
    public dtm2013
    public density_uncertainty_interface
    public density_uncertainty
    public P_ReadDTM12_interface
    public P_ReadDTM12

    interface dtm2013_interface
     subroutine dtm2013(day,f,fbar,akp,alti,hl,alat,xlon,tz,tinf,tp120,ro,d,wmm)
      real, intent(in) :: alat,alti,day,hl,xlon
      real, intent(out) :: ro,tinf,tz,tp120,wmm
      real, dimension(2), intent(in) :: f,fbar
      real, dimension(4), intent(in) :: akp
      real, intent(out) :: d(6)
     end subroutine
    end interface


    interface density_uncertainty_interface
     subroutine density_uncertainty(alt,lat,lst,flux,kp,unc)
      real, intent(in) :: lat,lst,kp
      real, intent(in) :: alt,flux
      real, intent(out) :: unc
     end subroutine
    end interface


    interface P_ReadDTM12_interface
     subroutine P_ReadDTM12()
     end subroutine
    end interface

end module dtm2013_interfaces
