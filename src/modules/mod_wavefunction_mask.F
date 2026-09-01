#if defined _GPU
#  include<dev_defs.h>
#else
#  include<dev_defs_nogpu.h>
#endif


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
      complex(SP), pointer :: wavefunction(:, :)
    end subroutine apply_WF_mask
  end interface


  interface  
    subroutine apply_WF_mask_gpu(wavefunction)
      use pars, only: SP
      complex(SP), pointer DEV_ATTR :: wavefunction(:, :)
    end subroutine apply_WF_mask_gpu
  end interface

  ! ...implementations in wfs_and_fft/wavefunction_mask_incl.F

end module wavefunction_mask_m


# include<dev_undefs.h>
