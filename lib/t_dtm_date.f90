!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                            !
!> @brief           Provide custom date type
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
!>   This module provides the dtm_date custom type. This is required to pass
!>   dates to the dtm_wrapper subroutine
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
module t_dtm_date
!> @brief Structure to hold the date input by the user
!> @details
!>This type allows the user to enter either a MJD2000 date or a day/month/year hh:mm:sec
!>date, without the need of passing a large number of parameters.
!>The user must set a type_flag and the required variables
!>Unused variables do not need to be initialized
type dtm_date
     integer :: type_flag !<1 for MJD2000 date, 2 for calendar date
     real*8  :: mjd2000 !< Date in MJD2000
     integer :: day !< Day
     integer :: month !< Month
     integer :: year !< Year
     integer :: hour !< Hour
     integer :: minute !< Minute
     real*8  :: second !< Seconds
end type
end module
