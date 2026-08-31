!
! License-Identifier: GPL
!
! Copyright (C) 2026 The Yambo Team
!
! Authors (see AUTHORS file for details): JK
!
#define WF_MASK
#define WF_SUB(x) CAT(x,_wavefunction_mask)
#include "scatter_Bamp_incl.F"
#undef WF_MASK
#define WF_SUB(x) x
#include "scatter_Bamp_incl.F"
#undef WF_SUB
