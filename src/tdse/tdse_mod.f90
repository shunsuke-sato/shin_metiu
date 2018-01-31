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
module tdse_mod
  use variables_for_tdse
  use model_parameters
  use math_parameters
  use io
  implicit none
  private

  public :: solve_tdse
  contains

    subroutine solve_tdse

      write(*,"(A)")"Start: Solving the exat time-dependent Schrodinger equation."
      call read_tdse_parameters

      write(*,"(A)")"Finish: Solving the exat time-dependent Schrodinger equation."
    end subroutine solve_tdse

    subroutine read_tdse_parameters
      integer :: iostat
      namelist/grid_parameters/ &
        lsize_ion, &
        lsize_elec, &
        nx_ion, &
        nx_elec


      rewind(fnum_input)
      read(fnum_input,nml=grid_parameters,iostat=iostat)



    end subroutine read_tdse_parameters


end module tdse_mod
