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
  use model_parameters
  use math_parameters
  use io
  implicit none
  private

  real(8) :: lsize_elec
  integer :: Nx_elec
  real(8) :: dx_elec

  real(8) :: lsize_ion
  integer :: Nx_ion
  real(8) :: dx_ion

  real(8),allocatable :: x_elec(:),x_ion(:)
  complex(8),allocatable :: zwfn(:,:)


  public :: solve_tdse
contains

  subroutine solve_tdse

    write(*,"(A)")"Start: Solving the exat time-dependent Schrodinger equation."
    call read_tdse_parameters
    
    call initialize_tdse
    
    write(*,"(A)")"Finish: Solving the exat time-dependent Schrodinger equation."
  end subroutine solve_tdse

  subroutine read_tdse_parameters
    integer :: iostat
    namelist/grid_parameters/ &
      lsize_ion, &
      lsize_elec, &
      nx_ion, &
      nx_elec
    
    write(*,"(A)")"Reading input parameters for TDSE calculation."
    
    rewind(fnum_input)
    read(fnum_input,nml=grid_parameters,iostat=iostat)
    
    write(fnum_inputlog,"(A)")"%grid_parameters"
    write(fnum_inputlog,"(A,2x,e26.16e3)")"lsize_elec =",lsize_elec
    write(fnum_inputlog,"(A,2x,I8)")      "nx_elec    =",nx_elec
    write(fnum_inputlog,"(A,2x,e26.16e3)")"lsize_ion =",lsize_ion
    write(fnum_inputlog,"(A,2x,I8)")      "nx_ion    =",nx_ion
    write(fnum_inputlog,"(A)")"/"
    
    
  end subroutine read_tdse_parameters
  
  subroutine initialize_tdse
    integer :: ix,iy
    real(8) :: xx,yy,ss
    
    write(*,"(A)")"Initialize the TDSE calculation."
    
    allocate(x_elec(0:nx_elec),x_ion(0:nx_ion))
    allocate(zwfn(0:nx_elec,0:nx_ion))
    
    dx_elec = lsize_elec/nx_elec
    do ix = 0,nx_elec
      x_elec(ix) = -0.5d0*lsize_elec + dx_elec*ix
    end do
    
    dx_ion = lsize_ion/nx_ion
    do ix = 0,nx_ion
      x_ion(ix) = -0.5d0*lsize_ion + dx_ion*ix
    end do

! temporal initial wave function
    do ix = 0, nx_ion
      xx = x_ion(ix)
      do iy = 0, nx_elec
        yy = x_elec(iy)
        zwfn(ix,iy) = exp(-xx**2)*exp(-yy**2)
      end do
    end do

    ss = sum(abs(zwfn)**2)*dx_ion*dx_elec
    zwfn = zwfn/sqrt(ss)
    
    
    
  end subroutine initialize_tdse
  
end module tdse_mod
