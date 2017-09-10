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
module io
  implicit none
  private

  character(16),parameter,public :: filename_input = './inp'
  integer,parameter,public       :: fnum_input     = 11

  character(16),parameter,public :: filename_inputlog = './inp.log'
  integer,parameter,public       :: fnum_inputlog     = 12



  public :: init_io, &
            end_io

contains
!-------------------------------------------------------------------------------
  subroutine init_io

    write(*,"(A)")"Initialize: I/O routines"
    open(fnum_input,file=filename_input)
    open(fnum_inputlog,file=filename_inputlog)

  end subroutine init_io
!-------------------------------------------------------------------------------
  subroutine end_io

    close(fnum_input)
    close(fnum_inputlog)
    write(*,"(A)")"Finalize: I/O routines"

  end subroutine end_io

end module io
