/**
 * This file is part of BSA Library.
 * Copyright (C) 2023  Michele Esposito Marzino 
 *
 * BSA Library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BSA Library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BSA Library.  If not, see <https://www.gnu.org/licenses/>.
 * */
#ifndef BSACL_BASE_H_
#define BSACL_BASE_H_

#ifndef REAL
# define REAL double
# define REAL_IS_DOUBLE
#endif
#define real REAL


#if (defined _DEBUG) && !(defined BSACLC_DEBUG)
# define BSACLC_DEBUG
#endif


#ifndef STRINGIFYMACRO_LITERAL
#  define STRINGIFYMACRO_LITERAL(X) #X
#  define xstr(s) STRINGIFYMACRO_LITERAL(s)
#  define STRINGIFYMACRO_VALUE(X) xstr(X)
#endif


#ifdef BSACL_USE_CONST_EXT_POINTERS
#  define __EXT_PTR_CONST const
#else
#  define __EXT_PTR_CONST
#endif


// Wind PSDs
#define BSACL_PSD_TYPE_UNKNOWN   0
#define BSACL_PSD_TYPE_VONKARMAN 1
#define BSACL_PSD_TYPE_KAIMAL    2
#define BSACL_PSD_TYPE_DAVENPORT 5



#ifdef BSACL_USE_CUDA__

#include <limits>

# define BSACL_PI (REAL)M_PI
# define ierr_t cudaError_t
# define BSACL_SUCCESS cudaSuccess
# define BSACL_DEVICE_FREE_MEM(X) cudaFree((void *)X)
# define BSACL_MEM
# define BSACL_INT int
# define BSACL_UINT unsigned int
# define BSACL_USHORT unsigned short
# define INT int
# define UINT unsigned int
# define USHORT unsigned short
# define BOOL BSACL_UINT
# define UINT_PTR_T unsigned int*
# define INT_PTR_T  int*
# define REAL_PTR_T REAL*
# define DEVICE __device__
# define PRIVATE
# define GLOBAL
# define LOCAL __shared__
# define CONSTANT
# define KERNEL __global__
# define EXTERN_C extern "C"
#ifdef REAL_IS_DOUBLE
# define POWR pow
#else
# define POWR powf
#endif
# define LOCAL_ID_X_DIM0 threadIdx.x
# define LOCAL_ID_Y_DIM1 threadIdx.y
# define LOCAL_ID_Z_DIM2 threadIdx.z
# define N_BLOCKS_WGROUPS_X_DIM0
# define N_BLOCKS_WGROUPS_Y_DIM1
# define N_BLOCKS_WGROUPS_Z_DIM2
# define BLOCK_ID_X_DIM0  blockIdx.x
# define BLOCK_ID_Y_DIM1  blockIdx.y
# define BLOCK_ID_Z_DIM2  blockIdx.z
# define GLOBAL_ID_X_DIM0 threadIdx.x + blockDim.x * BLOCK_ID_X_DIM0
# define GLOBAL_ID_Y_DIM1 threadIdx.y + blockDim.y * BLOCK_ID_Y_DIM1
# define GLOBAL_ID_Z_DIM2 threadIdx.z + blockDim.z * BLOCK_ID_Z_DIM2
# define LOCAL_WORKGROUP_BARRIER __syncthreads()

#else  // OpenCL

# define BSACL_PI      3.14159265358979323846
# define ierr_t cl_int
# define BSACL_SUCCESS CL_SUCCESS
# define BSACL_DEVICE_FREE_MEM(X) clReleaseMemObject(X)
# define BSACL_MEM cl_mem
# define BSACL_INT cl_int
# define BSACL_UINT cl_uint
# define BSACL_USHORT cl_uint
# define INT int
# define UINT uint
# define USHORT ushort
# define BOOL bool
# define UINT_PTR_T
# define INT_PTR_T
# define REAL_PTR_T
# define DEVICE
# define PRIVATE __private
# define GLOBAL __global
# define LOCAL __local
# define CONSTANT __constant
# define KERNEL __kernel
# define EXTERN_C
# define POWR powr
# define LOCAL_ID_X_DIM0 get_local_id(0)
# define LOCAL_ID_Y_DIM1 get_local_id(1)
# define LOCAL_ID_Z_DIM2 get_local_id(2)
# define N_BLOCKS_WGROUPS_X_DIM0 get_num_groups(0)
# define N_BLOCKS_WGROUPS_Y_DIM1 get_num_groups(1)
# define N_BLOCKS_WGROUPS_Z_DIM2 get_num_groups(2)
# define BLOCK_ID_X_DIM0  get_group_id(0)
# define BLOCK_ID_Y_DIM1  get_group_id(1)
# define BLOCK_ID_Z_DIM2  get_group_id(2)
# define GLOBAL_ID_X_DIM0 get_global_id(0)
# define GLOBAL_ID_Y_DIM1 get_global_id(1)
# define GLOBAL_ID_Z_DIM2 get_global_id(2)
# define LOCAL_WORKGROUP_BARRIER barrier(CLK_LOCAL_MEM_FENCE)

#endif


#ifdef REAL_IS_DOUBLE
# define REAL_MIN DBL_MIN
#else
# define REAL_MIN FLT_MIN
#endif


#endif // BSACL_BASE_H_
