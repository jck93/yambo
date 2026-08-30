#include <dev_defs.h>


module wavefunction_mask_m
  use pars, ONLY: SP
  real(SP), allocatable, target :: wf_mask(:)
  real(SP), allocatable, target DEV_ATTR :: wf_mask_d(:)


  interface  
    subroutine initialize_WF_mask()
    end subroutine initialize_WF_mask
  end interface


  interface  
    subroutine finalize_WF_mask()
    end subroutine finalize_WF_mask
  end interface


  interface  
    subroutine apply_WF_mask(wavefunction)
      use pars, only: SP
      complex(SP), pointer DEV_ATTR :: wavefunction(:, :)
    end subroutine apply_WF_mask
  end interface

  ! ...implementations in wfs_and_fft/wavefunction_mask.F

end module wavefunction_mask_m

