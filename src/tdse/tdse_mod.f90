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
  use global_variables
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
  real(8),allocatable :: v_pot(:,:)

  real(8),allocatable :: dwfn(:,:)
  complex(8),allocatable :: zwfn(:,:)

  real(8),allocatable :: dwfn_t(:,:),dhwfn_t(:,:)
  complex(8),allocatable :: zwfn_t(:,:),zhwfn_t(:,:)


  public :: solve_tdse
contains

  subroutine solve_tdse

    write(*,"(A)")"Start: Solving the exat time-dependent Schrodinger equation."
    call read_tdse_parameters
    
    call initialize_tdse

    select case(n_calc_mode)
      case(n_calc_mode_gs)
        call calc_ground_state_exact_schrodinger
      case(n_calc_mode_td)
        call calc_time_propagation_exact_schrodinger
      case default 
        write(*,"(A)")"Nothing computed in the exact Schrodinger equation mode."
      end select
    
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
    real(8) :: xx,yy,ss,rr
    real(8) :: func_tmp(0:nx_elec)
    
    write(*,"(A)")"Initialize the TDSE calculation."
    
    allocate(x_elec(0:nx_elec),x_ion(0:nx_ion))
    allocate(v_pot(0:nx_elec,0:nx_ion))
    allocate(dwfn(0:nx_elec,0:nx_ion))
    allocate(dwfn_t(0-2:nx_elec+2,0-2:nx_ion+2),dhwfn_t(0:nx_elec,0:nx_ion))
    dwfn_t = 0d0
      
    allocate(zwfn(0:nx_elec,0:nx_ion))

    allocate(zwfn_t(0-2:nx_elec+2,0-2:nx_ion+2),zhwfn_t(0:nx_elec,0:nx_ion))
    zwfn_t = 0d0
    
    dx_elec = lsize_elec/nx_elec
    do ix = 0,nx_elec
      x_elec(ix) = -0.5d0*lsize_elec + dx_elec*ix
    end do
    
    dx_ion = lsize_ion/nx_ion
    do ix = 0,nx_ion
      x_ion(ix) = -0.5d0*lsize_ion + dx_ion*ix
    end do

    do iy = 0, nx_ion
      yy = x_ion(iy)
      do ix = 0, nx_elec
        xx = x_elec(ix)

        v_pot(ix,iy) = 1d0/abs(yy-0.5d0*Ldist_m) + 1d0/abs(yy+0.5d0*Ldist_m)

        rr = (xx - 0.5d0*Ldist_m)/Rr_m
        v_pot(ix,iy) = v_pot(ix,iy) - erf_x(rr)/Rr_m

        rr = (xx + 0.5d0*Ldist_m)/Rl_m
        v_pot(ix,iy) = v_pot(ix,iy) - erf_x(rr)/Rl_m

        rr = (xx - yy)/Rf_m
        v_pot(ix,iy) = v_pot(ix,iy) - erf_x(rr)/Rf_m
          
      end do
    end do
    
    if(n_calc_mode == n_calc_mode_gs)then
      do iy = 0, nx_ion
        do ix = 0, nx_elec
          yy = x_ion(iy)
          xx = x_elec(ix)
          dwfn(iy,ix) = exp(-xx**2)*exp(-yy**2)
        end do
      end do

      ss = sum(abs(dwfn)**2)*dx_ion*dx_elec
      dwfn = dwfn/sqrt(ss)
    end if

    if(n_calc_mode == n_calc_mode_td)then
! temporal initial wave function
    open(500,file='init_wfn_elec.out', form='unformatted')
    read(500)ix
    if(ix /= nx_elec)stop 'Error: reading initial wnf_elec.'
    read(500)func_tmp(0:nx_elec)
    close(500)

    do iy = 0, nx_ion
      yy = x_ion(iy)+4d0
      zwfn(0:nx_elec,iy) = func_tmp(0:nx_elec)*exp(-0.5d0*(yy/sqrt(2.85d0))**2)
    end do

!    open(500,file='init_wfn_tot.out', form='unformatted')
!    read(500)dwfn
!    close(500)
!    zwfn = dwfn

    ss = sum(abs(zwfn)**2)*dx_ion*dx_elec
    zwfn = zwfn/sqrt(ss)

    end if
    
  end subroutine initialize_tdse

! zwfn_t => zhwfn_t
  subroutine calc_zhpsi_tdse
    real(8),parameter :: c0 = -5d0/2d0 &
      ,c1 =  4d0/3d0 &
      ,c2 = -1d0/12d0
    real(8) :: c0_e, c1_e, c2_e
    real(8) :: c0_i, c1_i, c2_i
    real(8) :: c0_ei
    integer :: ix,iy

    c0_e = -0.5d0*c0/dx_elec**2/mass_elec
    c1_e = -0.5d0*c1/dx_elec**2/mass_elec
    c2_e = -0.5d0*c2/dx_elec**2/mass_elec

    c0_i = -0.5d0*c0/dx_ion**2/mass_ion
    c1_i = -0.5d0*c1/dx_ion**2/mass_ion
    c2_i = -0.5d0*c2/dx_ion**2/mass_ion

    c0_ei = c0_e + c0_i

    do iy  = 0, nx_ion
      do ix = 0, nx_elec
        zhwfn_t(ix,iy) = (c0_ei + v_pot(ix,iy))*zwfn_t(ix,iy) &
          +c1_e*(zwfn_t(ix+1,iy) + zwfn_t(ix-1,iy)) &
          +c2_e*(zwfn_t(ix+2,iy) + zwfn_t(ix-2,iy)) &
          +c1_i*(zwfn_t(ix,iy+1) + zwfn_t(ix,iy-1)) &
          +c2_i*(zwfn_t(ix,iy+2) + zwfn_t(ix,iy-2))

      end do
    end do
    
  end subroutine calc_zhpsi_tdse

! dwfn_t => dhwfn_t
  subroutine calc_dhpsi_tdse
    real(8),parameter :: c0 = -5d0/2d0 &
      ,c1 =  4d0/3d0 &
      ,c2 = -1d0/12d0
    real(8) :: c0_e, c1_e, c2_e
    real(8) :: c0_i, c1_i, c2_i
    real(8) :: c0_ei
    integer :: ix,iy

    c0_e = -0.5d0*c0/dx_elec**2/mass_elec
    c1_e = -0.5d0*c1/dx_elec**2/mass_elec
    c2_e = -0.5d0*c2/dx_elec**2/mass_elec

    c0_i = -0.5d0*c0/dx_ion**2/mass_ion
    c1_i = -0.5d0*c1/dx_ion**2/mass_ion
    c2_i = -0.5d0*c2/dx_ion**2/mass_ion

    c0_ei = c0_e + c0_i

    do iy  = 0, nx_ion
      do ix = 0, nx_elec
        dhwfn_t(ix,iy) = (c0_ei + v_pot(ix,iy))*dwfn_t(ix,iy) &
          +c1_e*(dwfn_t(ix+1,iy) + dwfn_t(ix-1,iy)) &
          +c2_e*(dwfn_t(ix+2,iy) + dwfn_t(ix-2,iy)) &
          +c1_i*(dwfn_t(ix,iy+1) + dwfn_t(ix,iy-1)) &
          +c2_i*(dwfn_t(ix,iy+2) + dwfn_t(ix,iy-2))

      end do
    end do
    
  end subroutine calc_dhpsi_tdse

  subroutine dt_evolve_tdse(dt)
    real(8),intent(in) :: dt
    integer,parameter :: n_Taylor_order = 4
    integer :: iexp
    complex(8) :: zcoef

    zwfn_t(0:nx_elec,0:nx_ion) = zwfn(0:nx_elec,0:nx_ion)
    zcoef = 1d0

    do iexp = 1, n_Taylor_order
      zcoef = zcoef*(-zI*dt)/iexp
      
      call calc_zhpsi_tdse
      zwfn = zwfn + zcoef*zhwfn_t

      if(iexp == n_Taylor_order)exit
      zwfn_t(0:nx_elec,0:nx_ion) = zhwfn_t(0:nx_elec,0:nx_ion)

    end do
    

  end subroutine dt_evolve_tdse
  subroutine calc_time_propagation_exact_schrodinger
    integer :: nt, iter
    real(8) :: norm, Etot, entropy

    nt = aint(Tpropagation/time_step)+1
    write(*,"(A,2x,I8)")"# of time steps =",nt
    open(201,file='td_quantities.out')
    norm = sum(abs(zwfn)**2)*dx_elec*dx_ion
    zwfn_t(0:nx_elec,0:nx_ion) = zwfn(0:nx_elec,0:nx_ion)
    call calc_zhpsi_tdse
    call calc_von_neumann_entropy(entropy)
    Etot = sum(conjg(zwfn)*zhwfn_t)*dx_elec*dx_ion
    write(201,"(999e26.16e3)")0d0,norm,Etot,entropy
    write(*,"(999e26.16e3)")0d0,norm,Etot,entropy

    do iter = 1,nt
      write(*,"(A,2x,I7)")"iter=",iter
      call dt_evolve_tdse(time_step)

      if(mod(iter,100) == 0)then
        norm = sum(abs(zwfn)**2)*dx_elec*dx_ion
        zwfn_t(0:nx_elec,0:nx_ion) = zwfn(0:nx_elec,0:nx_ion)
        call calc_zhpsi_tdse
        call calc_von_neumann_entropy(entropy)
        Etot = sum(conjg(zwfn)*zhwfn_t)*dx_elec*dx_ion
        write(201,"(999e26.16e3)")time_step*iter,norm,Etot,entropy
        write(*,"(999e26.16e3)")time_step*iter,norm,Etot,entropy
      end if
    end do

    close(201)


  end subroutine calc_time_propagation_exact_schrodinger

  subroutine calc_ground_state_exact_schrodinger
    integer,parameter :: max_iter = 1000000
    real(8),parameter :: dt_imag = 0.01d0
    real(8) :: Etot, norm, res
    integer :: iter
    integer :: ix,iy
    real(8),allocatable :: rdm_ion(:,:)
!LAPACK
    integer :: nmax
    real(8),allocatable :: a_mat(:,:)
    integer :: lwork
    real(8),allocatable :: work_lp(:)
    real(8),allocatable :: rwork(:),w(:)
    integer :: info

    nmax = nx_ion+1
    lwork=6*nmax
    allocate(work_lp(lwork),rwork(3*nmax-2),w(nmax))

    allocate(a_mat(nmax,nmax))


    allocate(rdm_ion(0:nx_ion,0:nx_ion))

    do iter = 0,max_iter

      dwfn_t(0:nx_elec,0:nx_ion) = dwfn(0:nx_elec,0:nx_ion)
      call calc_dhpsi_tdse
      Etot = sum(dwfn*dhwfn_t)*dx_elec*dx_ion
      norm = sum(dwfn**2)*dx_elec*dx_ion
      Etot = Etot/norm
      res  = sum((dhwfn_t-Etot*dwfn)**2)*dx_elec*dx_ion
      res  = res/norm
      write(*,"(A,2x,I7,2x,2e26.16e3)")"iter, Etot, residual =",iter,Etot,res

      dwfn = dwfn - dt_imag * dhwfn_t
      norm = sum(dwfn**2)*dx_elec*dx_ion
      dwfn = dwfn/sqrt(norm)

      if(res < 1d-7)exit

    end do

    open(20,file='density_elec_gs.out')
    do ix = 0,nx_elec
      write(20,"(999e26.16e3)")x_elec(ix),sum(dwfn(ix,:)**2)*dx_ion
    end do
    close(20)

    open(20,file='density_ion_gs.out')
    do ix = 0,nx_ion
      write(20,"(999e26.16e3)")x_ion(ix),sum(dwfn(:,ix)**2)*dx_elec
    end do
    close(20)

    do ix = 0, nx_ion
      do iy = ix,nx_ion
        rdm_ion(ix,iy) = sum(dwfn(:,ix)*dwfn(:,iy))*dx_elec
        rdm_ion(iy,ix) = rdm_ion(ix,iy)
      end do
    end do

    norm = 0d0
    do ix = 0,nx_ion
      norm = norm + rdm_ion(ix,ix)
    end do
    norm = norm*dx_ion
    write(*,*)norm

    a_mat(1:nmax,1:nmax) = -rdm_ion(0:nx_ion,0:nx_ion)*dx_ion
    call dsyev('V', 'U', nmax, a_mat, nmax, w, work_lp, lwork, info)
    
    write(*,"(A)")"Occupation of natural ion-orbitals"
    write(*,*)sum(w)
    do ix = 1,20
      write(*,"(I7,2x,e26.16e3)")ix,-w(ix)
    end do

    open(20,file='natural_ion_orbital.out')
    do ix = 1,nmax
      write(20,"(999e26.16e3)")x_ion(ix-1),a_mat(ix,1:4)
    end do
    close(20)

    open(500,file='init_wfn_tot.out', form='unformatted')
    write(500)dwfn
    close(500)

  end subroutine calc_ground_state_exact_schrodinger

  subroutine calc_von_neumann_entropy(entropy_out)
    real(8),intent(out) :: entropy_out
    integer :: ix,iy
    real(8) :: norm
    complex(8),allocatable :: zrdm_ion(:,:)
!LAPACK
    complex(8),allocatable :: za_mat(:,:)
    integer :: nmax
    integer :: lwork
    complex(8),allocatable :: work_lp(:)
    real(8),allocatable :: rwork(:),w(:)
    integer :: info

    nmax = nx_ion+1
    lwork=6*nmax
    allocate(work_lp(lwork),rwork(3*nmax-2),w(nmax))

    allocate(za_mat(nmax,nmax))


    allocate(zrdm_ion(0:nx_ion,0:nx_ion))

    do ix = 0, nx_ion
      do iy = ix,nx_ion
        zrdm_ion(ix,iy) = sum(conjg(zwfn(:,ix))*zwfn(:,iy))*dx_elec
        zrdm_ion(iy,ix) = conjg(zrdm_ion(ix,iy))
      end do
    end do

    norm = 0d0
    do ix = 0,nx_ion
      norm = norm + zrdm_ion(ix,ix)
    end do
    norm = norm*dx_ion

    za_mat(1:nmax,1:nmax) = zrdm_ion(0:nx_ion,0:nx_ion)*dx_ion/norm

    Call zheev('V', 'U', nmax, za_mat, nmax, w, work_lp, lwork, rwork, info)

    entropy_out = 0d0
    do ix = 1,nmax
      if(w(ix) >0d0)entropy_out = entropy_out - w(ix)*log(w(ix))
    end do
    write(*,"(A,2x,e26.16e3)")"entropy=",entropy_out

  end subroutine calc_von_neumann_entropy
  
end module tdse_mod
