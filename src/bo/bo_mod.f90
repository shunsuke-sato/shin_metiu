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
  use model_parameters
  use math_parameters
  use io
  implicit none
  private


  real(8) :: lsize_elec
  integer :: Nx_elec
  real(8) :: Rion_bo
  integer :: nstate_bo

  real(8) :: dx_elec

  real(8),allocatable :: dwfn(:,:), vion(:), x_elec(:)

  public :: solve_bo
contains

  subroutine solve_bo

    write(*,"(A)")"Start: Solving the electronic static Schrodinger equation,"
    write(*,"(A)")"       under Born-Oppenheimer approximation."
    if(n_calc_mode /= n_calc_mode_gs)then
      error_messages(1) = "Error: The Born-Oppenheimer approximation should be used"
      error_messages(2) = "       with the ground state calculation mode."
      call error_finalize(2)
    end if

    call read_bo_parameters

    call initialize_bo

    call eigen_solver_bo
    
    write(*,"(A)")"Finish: Solving the electronic static Schrodinger equation,"
    write(*,"(A)")"       under Born-Oppenheimer approximation."
  end subroutine solve_bo
  
  subroutine read_bo_parameters
    integer :: iostat
    namelist/grid_parameters/ &
      lsize_elec, &
      nx_elec
    namelist/BO_parameters/ &
      Rion_bo, &
      nstate_bo
    
    write(*,"(A)")"Reading input parameters for BO calculation."

    
    rewind(fnum_input)
    read(fnum_input,nml=grid_parameters,iostat=iostat)

    rewind(fnum_input)
    read(fnum_input,nml=BO_parameters,iostat=iostat)

    write(fnum_inputlog,"(A)")"%grid_parameters"
    write(fnum_inputlog,"(A,2x,e26.16e3)")"lsize_elec =",lsize_elec
    write(fnum_inputlog,"(A,2x,I8)")      "nx_elec    =",nx_elec
    write(fnum_inputlog,"(A)")"/"

    write(fnum_inputlog,"(A)")"%BO_parameters"
    write(fnum_inputlog,"(A,2x,e26.16e3)")"Rion_bo     =",Rion_bo
    write(fnum_inputlog,"(A,2x,I8)")      "nx_state_bo =",nstate_bo
    write(fnum_inputlog,"(A)")"/"


    
  end subroutine read_bo_parameters

  subroutine initialize_bo
    integer :: ix
    real(8) :: xx

    write(*,"(A)")"Initialize the BO calculation."

    allocate(dwfn(0:nx_elec,nstate_bo), vion(0:nx_elec), x_elec(0:nx_elec))

    dx_elec = lsize_elec/nx_elec

    do ix = 0,nx_elec
      x_elec(ix) = -0.5d0*lsize_elec + dx_elec*ix
    end do

! ionic-potential
    vion = 0d0
    do ix = 0,nx_elec
      xx = (x_elec(ix) - 0.5d0*Ldist_m)/Rr_m
      vion(ix) = - erf_x(xx)/Rr_m

      xx = (x_elec(ix) + 0.5d0*Ldist_m)/Rl_m
      vion(ix) = vion(ix) - erf_x(xx)/Rl_m

      xx = (x_elec(ix) - Rion_bo)/Rf_m
      vion(ix) = vion(ix) - erf_x(xx)/Rf_m
    end do



  end subroutine initialize_bo

  subroutine eigen_solver_bo
    real(8),parameter :: c0 = -5d0/2d0 &
                        ,c1 =  4d0/3d0 &
                        ,c2 = -1d0/12d0
    real(8),allocatable :: a_mat(:,:)
    real(8) :: Eion_bo
!LAPACK
    integer :: nmax
    integer :: lwork
    real(8),allocatable :: work_lp(:)
    real(8),allocatable :: rwork(:),w(:)
    integer :: info
    integer :: ix, jx, ist

    nmax = nx_elec + 1
    lwork=6*nmax
    allocate(work_lp(lwork),rwork(3*nmax-2),w(nmax))

    allocate(a_mat(nmax,nmax))

    a_mat = 0d0
    do ix = 1, nmax
      do jx = 1,nmax
        if(abs(ix-jx)>2)then
        else if( abs(ix-jx) == 0 )then
          a_mat(ix,jx) = -0.5d0*c0/dx_elec**2
        else if( abs(ix-jx) == 1 )then
          a_mat(ix,jx) = -0.5d0*c1/dx_elec**2
        else if( abs(ix-jx) == 2 )then
          a_mat(ix,jx) = -0.5d0*c2/dx_elec**2
        end if
      end do
      a_mat(ix,ix) = a_mat(ix,ix) + vion(ix-1) 
    end do

    call dsyev('V', 'U', nx_elec, a_mat, nx_elec, w, work_lp, lwork, info)

    dwfn(0:nx_elec,1:nstate_bo) = a_mat(1:nmax,1:nstate_bo)
    Eion_bo = 1d0/abs(0.5d0*Ldist_m-Rion_bo) + 1d0/abs(0.5d0*Ldist_m+Rion_bo)
    write(*,"(A)")"States,  Electronic energy, Total energy"
    do ist = 1, nstate_bo
      write(*,"(I7,2x,2e26.16e3)")ist,w(ist),w(ist)+Eion_bo
    end do

    open(500,file='init_wfn_elec.out', form='unformatted')
    write(500)nx_elec
    write(500)dwfn(0:nx_elec,2)
    close(500)

    write(*,"(A)")"Solving eigenvalue problem for the BO calculation."
    
  end subroutine eigen_solver_bo

end module bo_mod
