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
module model_parameters
  use io
  implicit none
  private

  real(8),public :: Ldist_m
  real(8),public :: Rf_m, Rr_m, Rl_m
  real(8),public :: mass_elec, mass_ion

  public :: read_model_parameters

contains

  subroutine read_model_parameters
    integer :: iostat
    namelist/model_parameters/ &
      Ldist_m, &
      Rf_m, &
      Rr_m, &
      Rl_m, &
      mass_elec, &
      mass_ion


    rewind(fnum_input)
    read(fnum_input,nml=model_parameters,iostat=iostat)

    write(fnum_inputlog,"(A)")"%model_parameters"
    write(fnum_inputlog,"(A,2x,e26.16e3)")"Ldist_m =",Ldist_m
    write(fnum_inputlog,"(A,2x,e26.16e3)")"Rf_m =",Rf_m
    write(fnum_inputlog,"(A,2x,e26.16e3)")"Rr_m =",Rr_m
    write(fnum_inputlog,"(A,2x,e26.16e3)")"Rl_m =",Rl_m
    write(fnum_inputlog,"(A,2x,e26.16e3)")"mass_elec =",mass_elec
    write(fnum_inputlog,"(A,2x,e26.16e3)")"mass_ion =",mass_ion
    write(fnum_inputlog,"(A)")"/"

  end subroutine read_model_parameters
end module model_parameters
