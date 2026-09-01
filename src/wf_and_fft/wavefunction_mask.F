subroutine DEV_SUB(apply_WF_mask)(wavefunction)
  use pars, only: SP
  use FFT_m, only: fft_size
  use electrons, only: n_spinor
  use wavefunction_mask_m, only: WF_mask, WF_mask_d
  implicit none
  complex(SP), pointer DEV_ATTR :: wavefunction(:, :)
  real(SP), pointer DEV_ATTR :: WF_mask_p(:)
  real(SP) :: squared_norm, renormalization_factor
  logical :: wavefunction_is_spinor
  integer :: ir, grid_points
  integer, parameter :: spin_up = 1, spin_down = 2
  !
  !
  wavefunction_is_spinor = (n_spinor == 2)
  grid_points = fft_size
  squared_norm = 0.0_SP
  WF_mask_p => DEV_VAR(WF_mask)
  !
  ! First apply WF_mask to wavefunction and calculate its reduced norm...
  !
  !DEV_ACC_DEBUG data present(wavefunction, WF_mask_p)
  !DEV_ACC parallel loop reduction(+:squared_norm)
  !DEV_CUF kernel do(1) <<<*,*>>>
  !DEV_OMPGPU target map(tofrom: squared_norm) map(present, alloc: wavefunction, WF_mask_p)
  !DEV_OMPGPU teams loop reduction(+:squared_norm)
  !DEV_OMP parallel default(shared), private(ir), reduction(+:squared_norm)
  !DEV_OMP do
  do ir = 1, grid_points
    wavefunction(ir, spin_up) = wavefunction(ir, spin_up) * WF_mask_p(ir)
    squared_norm = squared_norm + abs(wavefunction(ir, spin_up)) ** 2
  enddo
  !DEV_OMPGPU end target
  !DEV_OMP end parallel    
  if (wavefunction_is_spinor) then
    !DEV_ACC_DEBUG data present(wavefunction, WF_mask_p)
    !DEV_ACC parallel loop reduction(+:squared_norm)
    !DEV_CUF kernel do(1) <<<*,*>>>
    !DEV_OMPGPU target map(tofrom:squared_norm) map(present, alloc: wavefunction, WF_mask_p)
    !DEV_OMPGPU teams loop reduction(+:squared_norm)
    !DEV_OMP parallel default(shared), private(ir), reduction(+:squared_norm)
    !DEV_OMP do
    do ir = 1, grid_points
      wavefunction(ir, spin_down) = wavefunction(ir, spin_down) * WF_mask_p(ir)
      squared_norm = squared_norm + abs(wavefunction(ir, spin_down)) ** 2
    enddo
    !DEV_OMPGPU end target
    !DEV_OMP end parallel    
  endif
  !
  ! ...then renormalize WF_masked wavefunction back to 1
  !
  renormalization_factor = 1.0_SP / sqrt(squared_norm)
  !DEV_ACC parallel loop copyin(renormalization_factor)
  !DEV_CUF kernel do(1) <<<*,*>>>
  !DEV_OMPGPU target map(to:renormalization_factor) map(present, alloc: wavefunction)
  !DEV_OMPGPU teams loop 
  !DEV_OMP parallel default(shared), private(ir)
  !DEV_OMP do 
  do ir = 1, grid_points
    wavefunction(ir, spin_up) = wavefunction(ir, spin_up) * renormalization_factor
  enddo
  !DEV_OMPGPU end target
  !DEV_OMP end parallel
  !DEV_ACC_DEBUG end data
  if (wavefunction_is_spinor) then
    !DEV_ACC parallel loop copyin(renormalization_factor)
    !DEV_CUF kernel do(1) <<<*,*>>>
    !DEV_OMPGPU target map(to: renormalization_factor) map(present, alloc: wavefunction)
    !DEV_OMPGPU teams loop 
    !DEV_OMP parallel default(shared), private(ir)
    !DEV_OMP do 
    do ir = 1, grid_points
      wavefunction(ir, spin_down) = wavefunction(ir, spin_down) * renormalization_factor
    enddo
    !DEV_OMPGPU end target
    !DEV_OMP end parallel
    !DEV_ACC_DEBUG end data
  endif
end subroutine DEV_SUB(apply_WF_mask)

