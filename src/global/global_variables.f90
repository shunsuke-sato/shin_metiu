!-------------------------------------------------------------------------------
!
!  Copyright 2017 Shunsuke A. Sato
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
module global_variables
  use io
  implicit none
  private

  integer,public :: n_calc_mode
  integer,public,parameter :: n_calc_mode_gs_td = 0
  integer,public,parameter :: n_calc_mode_gs    = 1
  integer,public,parameter :: n_calc_mode_td    = 2

  integer,public :: n_method
  integer,public,parameter :: n_method_exact = 0
  integer,public,parameter :: n_method_BO    = 1

  public :: read_global_variables

contains

  subroutine read_global_variables
    integer :: iostat
    character(256) :: calc_mode, method
    namelist/control/ &
      calc_mode, &
      method


    rewind(fnum_input)
    read(fnum_input,nml=control,iostat=iostat)

    write(fnum_inputlog,"(A)")"%control"
    write(fnum_inputlog,"(A,2x,A)")"calc_mode =",trim(calc_mode)
    write(fnum_inputlog,"(A,2x,A)")"method =",trim(method)
    write(fnum_inputlog,"(A)")"/"

    select case(calc_mode)
    case('gs_td','GS_TD')
       n_calc_mode = n_calc_mode_gs_td
    case('gs','GS')
       n_calc_mode = n_calc_mode_gs
    case('td','TD')
       n_calc_mode = n_calc_mode_td
    case default
       stop 'Error: invalid calc_mode'
    end select

    select case(method)
    case('exact')
       n_method = n_method_exact
    case('bo','BO')
       n_method = n_method_BO
    end select


  end subroutine read_global_variables
end module global_variables

