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
module bo_mod
  use global_variables
  use variables_for_bo
  use model_parameters
  use math_parameters
  use io
  implicit none
  private

  public :: solve_bo
contains

  subroutine solve_bo

    write(*,"(A)")"Start: Solving the electronic static Schrodinger equation,"
    write(*,"(A)")"       under Born-Oppenheimer approximation,"
    if(n_calc_mode /= n_calc_mode_gs)then
      error_messages(1) = "Error: The Born-Oppenheimer approximation should be used"
      error_messages(2) = "       with the ground state calculation mode."
      call error_finalize(2)
    end if
    call read_bo_parameters
    
    write(*,"(A)")"Finish: Solving the exat time-dependent Schrodinger equation."
  end subroutine solve_bo
  
  subroutine read_bo_parameters
    integer :: iostat
    namelist/grid_parameters/ &
      lsize_elec, &
      nx_elec
    
    
    rewind(fnum_input)
    read(fnum_input,nml=grid_parameters,iostat=iostat)
    

    
  end subroutine read_bo_parameters
  

end module bo_mod
