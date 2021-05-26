/*
  !---------------------------------------------------------------------!
  ! Created by Madu Manathunga on 04/29/2020                            !
  !                                                                     !                           
  ! Previous contributors: Yipu Miao                                    !
  !                                                                     !
  ! Copyright (C) 2020-2021 Merz lab                                    !
  ! Copyright (C) 2020-2021 Götz lab                                    !
  !                                                                     !
  ! This Source Code Form is subject to the terms of the Mozilla Public !
  ! License, v. 2.0. If a copy of the MPL was not distributed with this !
  ! file, You can obtain one at http://mozilla.org/MPL/2.0/.            !
  !_____________________________________________________________________!

  !---------------------------------------------------------------------!
  ! This source file contains device kernels pertaining to single       !
  ! precision ERI computation.                                          !
  !---------------------------------------------------------------------!
*/

#include "gpu.h"
#include <cuda.h>

#define SINGLE_PRECISION
#undef QUICKDouble
#define QUICKDouble float

static __constant__ int devTrans[TRANSDIM*TRANSDIM*TRANSDIM];
static __constant__ int Sumindex[10]={0,0,1,4,10,20,35,56,84,120};

#include "gpu_get2e_sp.h"
#include "gpu_get2e_subs_hrr_sp.h"
#include "int.h"


#define int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"
/*#include "gpu_get2e_subs_grad_sp.h"


//===================================

#undef int_spd
#define int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_grad_sp.h"

#undef int_spd
#undef int_spdf
#define int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_grad_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#define int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_grad_sp.h"
*/

#ifdef CUDA_SPDF
//===================================

#undef int_spd
#define int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#define int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"


#undef int_spd
#undef int_spdf
#undef int_spdf2
#define int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"


#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#define int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#define int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"



#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#define int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#define int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"



#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#define int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"


#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#define int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#define int_spdf10
#include "gpu_get2e_subs_sp.h"
#endif

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10

//Include the kernels for open shell eri calculations
#define OSHELL
#define int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#undef new_quick_2_gpu_get2e_subs_h
#include "gpu_get2e_subs_sp.h"
/*#include "gpu_get2e_subs_grad_sp.h"

//===================================

#undef int_spd
#define int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_grad_sp.h"

#undef int_spd
#undef int_spdf
#define int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_grad_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#define int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_grad_sp.h"
*/

#ifdef CUDA_SPDF
//===================================

#undef int_spd
#define int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#define int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#define int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#define int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#define int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#define int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#define int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#define int_spdf8
#undef int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#define int_spdf9
#undef int_spdf10
#include "gpu_get2e_subs_sp.h"

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#define int_spdf10
#include "gpu_get2e_subs_sp.h"
#endif

#undef int_spd
#undef int_spdf
#undef int_spdf2
#undef int_spdf3
#undef int_spdf4
#undef int_spdf5
#undef int_spdf6
#undef int_spdf7
#undef int_spdf8
#undef int_spdf9
#undef int_spdf10

#undef OSHELL
#undef SINGLE_PRECISION
#undef QUICKDouble
#define QUICKDouble double
