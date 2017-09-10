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
module control_parameters
  use io
  implicit none
  private

  character(8),public :: method
  
  public :: read_control_parameters

contains
  subroutine read_control_parameters
    integer :: iostat
    namelist/control_parameters/ &
      method

    rewind(fnum_input)
    read(fnum_input,nml=control_parameters,iostat=iostat)
    write(fnum_inputlog,"(A)")"%control_parameters"
    write(fnum_inputlog,"(A,2x,A)")"method =",trim(method)
    write(fnum_inputlog,"(A)")"/"
    
  end subroutine read_control_parameters


end module control_parameters
