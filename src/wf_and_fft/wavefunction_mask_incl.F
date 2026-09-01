!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details): JK
!
!


subroutine initialize_WF_mask()
  use gpu_m, only: have_gpu
  use wavefunction_mask_m, only: WF_mask, WF_mask_d
  call allocate_WF_mask()
  call fill_WF_mask()

contains

  subroutine allocate_WF_mask()
    use FFT_m, only: FFT_size
#   include <y_memory.h>
    YAMBO_ALLOC(WF_mask, (FFT_size))
    if (have_gpu) then
      YAMBO_ALLOC_GPU(WF_mask_d, (FFT_size))
    endif
  end subroutine allocate_WF_mask


  subroutine fill_WF_mask()
    use pars, only: SP
    use FFT_m, only: FFT_dim
    use devxlib, only: devxlib_memcpy_h2d
    use D_lattice, only: a, atom_pos, atom_mass
    implicit none
    integer :: n(3), i, i1, i2, i3
    integer :: is, ia
    real(SP) :: r(3)
    !
    !
    n = FFT_dim
    WF_mask = 0.0_SP
    i = 1
    do i3 = 0, n(3)-1
      do i2 = 0, n(2)-1
        do i1 = 0, n(1)-1
          r = position_vector(i1, i2, i3)
          if (is_inside_cavity(r)) WF_mask(i) = 1.0_SP
          i = i + 1
        end do
      end do
    end do
    if (have_gpu) call devxlib_memcpy_h2d(WF_mask_d, WF_mask)
  end subroutine fill_WF_mask


  function position_vector(i1, i2, i3) result(r)
    use pars, only: SP
    use FFT_m, only: FFT_dim
    use D_lattice, only: a
    integer :: i1, i2, i3
    integer :: n(3)
    real(SP) :: r(3)
    n = FFT_dim
    r = a(1, :)*i1/n(1) + a(2, :)*i2/n(2) + a(3, :)*i3/n(3)
  end function position_vector


  logical function is_inside_cavity(r)
    use pars, only: SP
    use D_lattice, only: atom_pos, n_atoms_species, n_atomic_species, Z_species
    implicit none
    real(SP) :: r(3)
    real(SP) :: dr(3)
    integer :: species, atom, Z
    real(SP) :: square_sphere_radius
    real(SP), parameter :: scaling_factor = 1.1
    !
    !
    is_inside_cavity = .false.
    do species = 1, n_atomic_species
      Z = Z_species(species)
      square_sphere_radius = (vdW_radius(Z) * scaling_factor) ** 2
      do atom = 1, n_atoms_species(species)
        dr = r - atom_pos(:, atom, species) 
        if (dot_product(dr, dr) < square_sphere_radius) then
          is_inside_cavity = .true.
          return
        end if
      end do
    end do
  end function is_inside_cavity


  real(SP) function vdW_radius(Z)
    use pars, only: SP
    use UNITS, only: BO2ANG
    integer :: Z
    real(SP), parameter :: vdW_radius_angstrom(54) = [ &
      1.001, 1.012, 0.825, 1.408, 1.485, 2.000, 1.583, 1.500, &
      1.287, 1.243, 1.144, 1.364, 1.639, 1.716, 1.705, 1.683, &
      1.639, 1.595, 1.485, 1.474, 1.562, 1.562, 1.562, 1.562, &
      1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.650, 1.727, &
      1.760, 1.771, 1.749, 1.727, 1.628, 1.606, 1.639, 1.639, &
      1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, &
      1.672, 1.804, 1.881, 1.892, 1.892, 1.881 ]
    vdW_radius = vdW_radius_angstrom(Z) / BO2ANG
  end function vdW_radius
end subroutine initialize_WF_mask


subroutine finalize_WF_mask()
  use gpu_m, only: have_gpu
  use wavefunction_mask_m, only: WF_mask, WF_mask_d
# include <y_memory.h>
  YAMBO_FREE(WF_mask)
  if (have_gpu) then
    YAMBO_FREE_GPU(WF_mask_d)
  endif
end subroutine finalize_WF_mask


! apply_WF_mask subroutine has a CPU and a GPU version 
!
! CPU
!
#include<dev_defs_nogpu.h>
#undef _GPU_LOC
#include"wavefunction_mask.F"
#include<dev_undefs.h>
!
! GPU; this includes a host-fallback when GPUs are not available
!
#if defined _GPU
#  include<dev_defs.h>
#  define _GPU_LOC
#  include"wavefunction_mask.F"
#else
#  include<dev_defs_nogpu.h>
#  undef _GPU_LOC
#  undef DEV_SUB
#  define DEV_SUB(x) CAT(x,_gpu)
#  include"wavefunction_mask.F"
#  include<dev_undefs.h>
#endif

