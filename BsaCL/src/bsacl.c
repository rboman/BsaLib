
//-------------------------------------
// General  MACROs
//-------------------------------------
#include "_base.h"
#ifndef BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP
#  define BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP 64
#endif

#ifndef BSACL_OPT_N_WORK_GROUPS
#  define BSACL_OPT_N_WORK_GROUPS 10
#endif

#ifndef MIN_N_OF_WORK_GROUPS
#  define MIN_N_OF_WORK_GROUPS 2
#endif

#ifndef BSACL_MAX_GPU_COUNT
#  define BSACL_MAX_GPU_COUNT 8
#endif

#ifndef BSACL_KERNEL_ID
#  define BSACL_KERNEL_ID 3
#endif

#ifdef BSACL_USE_CUDA__
# define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef BASE_DIRECTORY
# define BASE_DIRECTORY "./"  // By default, current binary directory.
#endif


#include "bsacl.h"
#ifdef BSACL_USE_CUDA__
# ifdef BSACL_PASS_PARAMS_BY_MACRO
#  undef BSACL_PASS_PARAMS_BY_MACRO  // in CUDA, we pass always by kernel param!
# endif

# define BSACL_CONV_PULSATION

// # define __CUDACC_RTC__
# pragma message("   --- [NOTE]:  Using  CUDA  in place of  OpenCL  specification!")
# include "bsacl.cl"
# ifdef NULL // CUDA does not like (void *) ptr assignment
#   undef  NULL
#   define NULL 0
# endif
#else // using OpenCL
# define CL_USE_DEPRECATED_OPENCL_1_2_APIS
# pragma message("   --- [NOTE]:  Using  OpenCL  specification!")
# define __CL_ENABLE_EXCEPTIONS
# include <CL/cl.h>
# define BSACL_GENERATE_PTX_BINARY
#endif





//-------------------------------------
// GLOBAL OpenCL VARIABLES
//-------------------------------------
unsigned int bsacl_init_called__ = 0;

#ifndef BSACL_USE_CUDA__

cl_uint         n_platforms__  = 0;
cl_platform_id *platform_all__ = NULL; // array of all available platforms
cl_platform_id  platform__     = NULL; // final selected platform

#  ifdef BSACLC_DEBUG
cl_device_type device_type__ = CL_DEVICE_TYPE_CPU;
#  else
cl_device_type device_type__ = CL_DEVICE_TYPE_DEFAULT;
#  endif
cl_uint        n_devices__   = 0;       // n. devices for finally selected platform
cl_device_id   *devices__    = NULL;    // array of all available devices for selected platform
unsigned short i_def_dev_id_ = 0U;

/* Max memory allocation (bytes) allowed for one single buffer/image */
cl_ulong *dev_max_mem_alloc_size__  = NULL;
cl_ulong *dev_glob_mem_size__       = NULL;
cl_ulong *dev_glob_mem_cache_size__ = NULL;
cl_ulong *dev_local_mem_size__      = NULL;
cl_uint  *dev_n_compute_units__     = NULL;
size_t   *dev_max_work_group_size__ = NULL;
cl_uint  *dev_max_work_items_n_dims__ = NULL;
size_t  **dev_max_work_item_sizes__   = NULL;

/* Context */
cl_context context__ = NULL;

/* Array of devices' command queues */
cl_command_queue *cqueues__ = NULL;


cl_program program__;
char     **clSourceStrings__ = NULL;
char     **evalFctStrings__  = NULL;
cl_kernel  kernel_bfm__;
char       prog_compile_opts__[BUFSIZ] = {' '};

#else

dim3 cuda_nblocks__, cuda_nthreads_perblock__;

#endif // BSACL_USE_CUDA__


BSACL_UINT n_work_dims__;
size_t     global_dims_ie_NTWI__[3];
size_t     local_dims_ie_WIpWG__[3];
size_t     n_work_groups__[3];


// BSACL memory buffers
BSACL_MEM  UINT_PTR_T  d_tc__           = NULL;
BSACL_MEM  UINT_PTR_T  d_nodes_load__   = NULL;
BSACL_MEM  REAL_PTR_T  d_phiTc__        = NULL;
BSACL_MEM  REAL_PTR_T  d_nod_corr__     = NULL;
#if (BSACL_KERNEL_ID==1)
BSACL_MEM  REAL_PTR_T  d_S_uvw_fi__     = NULL;
BSACL_MEM  REAL_PTR_T  d_S_uvw_fj__     = NULL;
BSACL_MEM  REAL_PTR_T  d_S_uvw_fiPfj__  = NULL;
BSACL_MEM  REAL_PTR_T  d_fi_abs__       = NULL;
BSACL_MEM  REAL_PTR_T  d_fj_abs__       = NULL;
BSACL_MEM  REAL_PTR_T  d_fiPfj_abs__    = NULL;
#elif (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
BSACL_MEM  REAL_PTR_T  d_wind_nodal_vel__   = NULL;
BSACL_MEM  REAL_PTR_T  d_wind_turb_scales__ = NULL;
BSACL_MEM  REAL_PTR_T  d_wind_turb_std__    = NULL;
BSACL_MEM  INT_PTR_T   d_wind_nodal_windz__ = NULL;
#  if (BSACL_KERNEL_ID==3)
BSACL_MEM  REAL_PTR_T  d_fi__ = NULL;
BSACL_MEM  REAL_PTR_T  d_fj__ = NULL;
#  endif
#endif
BSACL_MEM  REAL_PTR_T  d_m3mf__ = NULL;



// ------------------------------------------------------------------
// EXTERNAL DATA (to be acquired, i.e. no internal memory allocation)
// ------------------------------------------------------------------

extdata_t extdata__;
#if (BSACL_KERNEL_ID==1)
unsigned int BASE_PSD_ALLOC_SIZE__ = 0U;
evalFct_t evalFunc_ptr__ = NULL;
#endif
unsigned int nfi__, nfj__;
REAL *fi__ = NULL;
REAL *fj__ = NULL;
REAL *S_uvw__       = NULL;
REAL *S_uvw_fiPfj__ = NULL;

unsigned short has_halted__  = 0U;
unsigned short has_cleaned__ = 0U;


#ifdef __cplusplus
extern "C" {
#endif


// -----------------------------------
//  LOCAL helper functions
// -----------------------------------


ierr_t releaseDeviceMem_() {
   unsigned ierr_;

   ierr_  = BSACL_DEVICE_FREE_MEM(d_m3mf__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_nodes_load__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_nod_corr__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_phiTc__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_tc__);

#if (BSACL_KERNEL_ID==1)   
   ierr_ != BSACL_DEVICE_FREE_MEM(d_fi_abs__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_fj_abs__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_fiPfj_abs__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_S_uvw_fi__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_S_uvw_fiPfj__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_S_uvw_fj__);
#elif (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_nodal_vel__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_nodal_windz__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_turb_scales__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_turb_std__);
# if (BSACL_KERNEL_ID==3)
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_fi__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_fj__);
# endif
#endif // (BSACL_KERNEL_ID==2)
   return (ierr_t)ierr_;
}

/**
 * @brief releases local memory
 * 
 */
void freeMem_(void) {

   if (has_cleaned__) return;

#ifndef BSACL_USE_CUDA__

   if (clSourceStrings__ != NULL) {
      free(*clSourceStrings__);
      clSourceStrings__ = NULL;
   }

   if (platform_all__ != NULL) free(platform_all__);
   platform_all__ = NULL;
   platform__     = NULL;
   
   if (devices__ != NULL) free(devices__);
   devices__ = NULL;

   if (dev_max_mem_alloc_size__ != NULL) free(dev_max_mem_alloc_size__);
   dev_max_mem_alloc_size__ = NULL;
   if (dev_glob_mem_size__ != NULL)      free(dev_glob_mem_size__);
   dev_glob_mem_size__ = NULL;
   if (dev_glob_mem_cache_size__ != NULL) free(dev_glob_mem_cache_size__);
   dev_glob_mem_cache_size__ = NULL;
   if (dev_local_mem_size__ != NULL)      free(dev_local_mem_size__);
   dev_local_mem_size__ = NULL;
   if (dev_n_compute_units__ != NULL)     free(dev_n_compute_units__);
   dev_n_compute_units__ = NULL;
   if (dev_max_work_group_size__ != NULL)   free(dev_max_work_group_size__);
   dev_max_work_group_size__ = NULL;
   if (dev_max_work_items_n_dims__ != NULL) free(dev_max_work_items_n_dims__);
   dev_max_work_items_n_dims__ = NULL;

   if (dev_max_work_item_sizes__ != NULL) {
      for (unsigned int i = 0; i < n_devices__; ++i) {
         if (dev_max_work_item_sizes__[i] != NULL) free(dev_max_work_item_sizes__[i]);
      }
      free(dev_max_work_item_sizes__);
   }

   if (cqueues__ != NULL) free(cqueues__);
   cqueues__ = NULL;
#endif

   freeExtData(&extdata__);

#if (BSACL_KERNEL_ID==1)
   evalFunc_ptr__ = NULL;
#endif

   fi__    = NULL;
   fj__    = NULL;
   S_uvw__ = NULL;

   has_cleaned__       = 1U;
   bsacl_init_called__ = 0U;
}



void printMsg_(const char *const msgtype, const char *const msg) {
   printf("\n%s%s.", msgtype, msg);
}
void printMsgWithIerr_(const char *const msgtype, const char *const msg, const int ierr) {
   printMsg_(msgtype, msg);
   if (ierr != 0) printf(" Aborting  (%d)", ierr);
}
void abortInternal_(const int ierr, const char *const emsg) {
   if (emsg != NULL) printMsgWithIerr_(ERRR_MSG, emsg, ierr);
   freeMem_();
   has_halted__ = 1U;
}


#define ABORT_INTERNAL_RETURN_VOID(X, Y)  {abortInternal_(X, Y); return;}
#define ABORT_INTERNAL_GOTO_RET_(X, Y, Z) {abortInternal_(X, Y); goto Z;}
#define ABORT_INTERNAL_RETURN_IERR(X, Y)  {abortInternal_(X, Y); return (X);}





#ifndef BSACL_USE_CUDA__
/**
 * @brief Get the all OpenCL Platforms available
 * 
 * @param npltfms pointer to the total n. of platforms found
 * @return BSACL_INT : error code
 */
BSACL_INT getNumPlatforms_(cl_uint *const npltfms) {
   BSACL_INT ierr_;

#ifdef BSACLC_DEBUG
   printMsg_(INFO_MSG, "Init querying for n. of available platforms...");
#endif

   ierr_ = clGetPlatformIDs(0, NULL, npltfms);
   if (*npltfms <= 0) {
      printMsg_(ERRR_MSG, "N<=0. OpenCL platform could not be correctly found.");
      return BSACL_PLATFORM_FIND_ERROR_;
   }
#ifdef BSACLC_DEBUG
   printf("\n%sFound  %d  platforms.", INFO_MSG, *npltfms);
#endif
   return ierr_;
}




BSACL_INT getAllPlatforms_() {
   BSACL_INT ierr_;

   ierr_ = getNumPlatforms_(&n_platforms__);
   if (ierr_ != BSACL_SUCCESS) {
      printf("\n%sCannot inquire n. of available CL platforms. Aborting.", ERRR_MSG);
      return ierr_;
   }

   platform_all__ = (cl_platform_id *) malloc(sizeof(cl_platform_id) * n_platforms__);
   ierr_          = clGetPlatformIDs(n_platforms__, platform_all__, NULL);
   if (ierr_ != BSACL_SUCCESS) {
      printf("\n%sCannot acquire  %d  available platforms. Aborting.", ERRR_MSG, n_platforms__);
      return ierr_;
   }
#ifdef BSACLC_DEBUG
   printf("\n%s%d  platforms acquired -- ok.", INFO_MSG, n_platforms__);
#endif
   return ierr_;
}



/**
 * @brief Select a platform by given platform info
 * 
 * @param selected_platform pointer to the actual platform criteria
 * @param pinfo platform info to query
 * @param pinfo_val value of info to query (if needed)
 * @param pinfo_nbytes n bytes of pinfo value (if any)
 * @return BSACL_INT 
 */
BSACL_INT getCLPlatformFromInfo_(
      cl_platform_id *selected_platform, cl_platform_info pinfo,
      const char * const pinfo_val,      size_t pinfo_nbytes)
{
   BSACL_INT ierr_;
   size_t param_size_;
   char buf_[BUFSIZ];

   *selected_platform = NULL;  // default initialise final platform to null

   if (platform_all__ == NULL) {
      ierr_ = getAllPlatforms_();
      if (platform_all__ == NULL) ierr_ = (BSACL_INT)-1;
      if (ierr_ != BSACL_SUCCESS) return ierr_;
   }

   if (pinfo_val != NULL) { // Search for info requested (all char[] !!)

      for (unsigned int i = 0; i < n_platforms__; i++) {

         ierr_ = clGetPlatformInfo(
            platform_all__[i], pinfo, BUFSIZ, &buf_, &param_size_);
         if (ierr_ == BSACL_SUCCESS) {
            if (strstr(buf_, pinfo_val) != NULL) {
               printf("\n%sPlatform with info  %s  found.", NOTE_MSG, pinfo_val);
               *selected_platform = platform_all__[i];
               break;
            }
         }
      }
   }

   // If info NOT found, or val not provided
   // Acquire first by default
   if (*selected_platform == NULL) {

      *selected_platform = platform_all__[0];

      if (*selected_platform == NULL) {
         printMsg_(ERRR_MSG, "Cannot even acquire first default CL platform. Check.");
         return BSACL_DEF_PLATFORM_ACQUIRE_ERROR_;
      }
      ierr_ = clGetPlatformInfo(
         *selected_platform, CL_PLATFORM_NAME, 1024, &buf_, &param_size_);
      
      buf_[param_size_] = '\0'; // BUG: make sure string is NULL terminated
      printf("\n%sPlatform with info   %s   NOT FOUND.", WARN_MSG, pinfo_val);
      printf("\n%sSelecting first by default (name: %s)", CONT_MSG, buf_);
   }
   printf("\n");
   return BSACL_SUCCESS;
};



/**
 * @brief Get the Selected Devices critical info needed for optimal setting.
 * 
 * @param  devices  list of selected available devices
 * @return (cl_int) Error code 
 */
void getCLDevicesInfo_() {
   size_t sz_, info_len_;

   // Allocate memory
   sz_ = sizeof(cl_ulong) * n_devices__;
   dev_max_mem_alloc_size__  = (cl_ulong *) malloc(sz_);
   dev_glob_mem_size__       = (cl_ulong *) malloc(sz_);
   dev_glob_mem_cache_size__ = (cl_ulong *) malloc(sz_);
   dev_local_mem_size__      = (cl_ulong *) malloc(sz_);

   sz_ = sizeof(cl_uint)  * n_devices__;
   dev_n_compute_units__       = (cl_uint *) malloc(sz_);
   dev_max_work_items_n_dims__ = (cl_uint *) malloc(sz_);
   sz_ = sizeof(size_t)  * n_devices__;
   dev_max_work_group_size__   = (size_t *) malloc(sz_);

   sz_ = sizeof(size_t*) * n_devices__;
   dev_max_work_item_sizes__ = (size_t **) malloc(sz_);

   for (unsigned int i = 0; i < n_devices__; ++i) {

#ifdef BSACLC_DEBUG
      printf("\n\n%sCritical info for device   #%d\n", INFO_MSG, i+1);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), (void*)&(dev_glob_mem_size__[i]), 
         &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - glob mem size       =   %llu\n", CONT_MSG, dev_glob_mem_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(cl_ulong), 
         (void*)&dev_glob_mem_cache_size__[i], &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - glob mem cache size =   %llu\n", CONT_MSG, dev_glob_mem_cache_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), 
         (void*)&dev_local_mem_size__[i], &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - local mem size      =   %llu\n", CONT_MSG, dev_local_mem_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), 
         (void*)&dev_max_mem_alloc_size__[i], &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - max mem alloc size  =   %llu\n", CONT_MSG, dev_max_mem_alloc_size__[i]);
#endif


      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), 
         (void*)&dev_n_compute_units__[i], &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - max compute units   =   %u\n", CONT_MSG, dev_n_compute_units__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), 
         (void*)&dev_max_work_group_size__[i], &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - max work-group size =   %llu\n", CONT_MSG, dev_max_work_group_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), 
         (void*)&dev_max_work_items_n_dims__[i], &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - max work-item dims  =   %u\n", CONT_MSG, dev_max_work_items_n_dims__[i]);
#endif

      sz_ = sizeof(size_t) * dev_max_work_items_n_dims__[i];
      dev_max_work_item_sizes__[i] = (size_t *) malloc(sz_);
      memset(dev_max_work_item_sizes__[i], (size_t)0, sz_);
      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_WORK_ITEM_SIZES, 
         sizeof(size_t) * dev_max_work_items_n_dims__[i], 
         (void*)dev_max_work_item_sizes__[i], &info_len_);
#ifdef BSACLC_DEBUG
      printf("%s - max work-item sizes =   ", CONT_MSG);
      for (unsigned int j = 0; j < dev_max_work_items_n_dims__[i] - 1; ++j)
         printf("%llu - ", dev_max_work_item_sizes__[i][j]);
      printf("%llu\n", dev_max_work_item_sizes__[i][dev_max_work_items_n_dims__[i] - 1]);
#endif

   }
}





/**
 * @brief Get the Devices By Type object
 * 
 * @param platform Chosen platform
 * @param device_type Type of desired devices to query for.
 * @param devices  Stores the actual devices (of type) found for platform
 * @return (cl_int) Number of found devices (of specified type)
 */
cl_int getCLDevicesByType_(
      cl_platform_id *platform, cl_device_type device_type, cl_device_id **devices) {

   cl_int ierr_ = 0;
   cl_uint num_devices_of_type_ = 0;

#ifdef BSACLC_DEBUG
   printMsg_(INFO_MSG, "Init querying devices by type..");
#endif

   // Get all available devices of type, for given platform
   ierr_ = clGetDeviceIDs(*platform, device_type, 0, NULL, &num_devices_of_type_);
   if (ierr_ != BSACL_SUCCESS) {
      printf(
         "\n\n%sCannot query requested devices type for provided platform. Aborting.", ERRR_MSG);
      return -1;
   }
   *devices = (cl_device_id *)malloc(sizeof(cl_device_id) * num_devices_of_type_);
   if (*devices == NULL) {
      printf("\n\n%sFailed to allocate memory for type of devices. Aborting.", ERRR_MSG);
      return -2;
   }

   // get actual devices by given type
   ierr_ = clGetDeviceIDs(*platform, device_type, num_devices_of_type_, *devices, NULL);
   if (ierr_ != BSACL_SUCCESS || devices == NULL) {
      printf("\n\n%sFailed to acquire type of devices for platform. Aborting (%d).", ERRR_MSG, ierr_);
      return -3;
   }
#ifdef BSACLC_DEBUG
   printf("\n%sFound  %d  device(s) of given type.", INFO_MSG, num_devices_of_type_);
#endif
   return num_devices_of_type_;
}
#endif // #ifndef BSACL_USE_CUDA__




void initCreateDeviceBuffers_() {
   ierr_t ierr_;
   size_t sz_;

   if (extdata__.tc__ == NULL)
      ABORT_INTERNAL_RETURN_VOID(-83, "List of turb. components was not acquired. Aborting.");
   sz_ = extdata__.NTC__ * sizeof(unsigned int);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_tc__, sz_);
#else
   d_tc__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.tc__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_tc__==NULL)
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create  d_tc__  device buffer");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_tc__, (void *)extdata__.tc__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_tc__.");
#endif


   if (extdata__.nodes_load__ == NULL)
      ABORT_INTERNAL_RETURN_VOID(-83, "List of loaded nodes was not acquired. Aborting.");
   sz_ = extdata__.NNODES_LOAD__ * sizeof(unsigned int);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_nodes_load__, sz_);
#else
   d_nodes_load__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.nodes_load__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_nodes_load__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create  d_nodes_load__  device buffer");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_nodes_load__, (void *)extdata__.nodes_load__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_nodes_load__.");
#endif



   if (extdata__.phi_T_c__ == NULL)
      ABORT_INTERNAL_RETURN_VOID(-83, "phiTc was not acquired. Aborting.");
   sz_ = extdata__.NMODES_EFF__ * extdata__.NNODES_LOAD__ * extdata__.NDEGW__ * sizeof(REAL);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_phiTc__, sz_);
#else
   d_phiTc__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.phi_T_c__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_phiTc__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create phi_T_c__ on device");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_phiTc__, (void *)extdata__.phi_T_c__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_phiTc__.");
#endif

   if (extdata__.nod_corr__ == NULL)
      ABORT_INTERNAL_RETURN_VOID(-83, "Nodal spatial correlation was not acquired. Aborting.");
   sz_ = extdata__.NNOD_CORR__ * extdata__.NTC__ * sizeof(REAL);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_nod_corr__, sz_);
#else
   d_nod_corr__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.nod_corr__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_nod_corr__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_nod_corr__ on device");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_nod_corr__, (void *)extdata__.nod_corr__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_nod_corr__.");
#endif


#if (BSACL_KERNEL_ID==1)

   // BUG: check if not using CL_MEM_COPY_HOST_PTR ?!
   // NOTE: here do not use any host ptr, enqueue write manually
   if (S_uvw__ == NULL)
      ABORT_INTERNAL_RETURN_VOID(-83, "Computed PSD values were not acquired. Aborting.");
   sz_ = BASE_PSD_ALLOC_SIZE__ * sizeof(REAL);
   // d_S_uvw_fi__ = clCreateBuffer(
   //    context__, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, 
   //    sz_, (void *)S_uvw__, &ierr_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_S_uvw_fi__, sz_);
#else
   d_S_uvw_fi__ = clCreateBuffer(context__, CL_MEM_READ_ONLY, sz_, NULL, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_S_uvw_fi__ device buffer");

#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_S_uvw_fj__, sz_);
#else
   d_S_uvw_fj__ = clCreateBuffer(context__, CL_MEM_READ_ONLY, sz_, NULL, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_S_uvw_fj__ device buffer");

   S_uvw_fiPfj__ = (REAL *)malloc(sz_);
   if (S_uvw_fiPfj__ == NULL)
      ABORT_INTERNAL_RETURN_VOID(-432, "Cannot allocate  HOST  memory for S_uvw_fiPfj__.");
// #ifdef BSACLC_DEBUG
//    else
//       printf("\n%sAllocated  S_uvw_fiPfj__  with size = %d", DEBG_MSG, sz_ / sizeof(REAL));
// #endif
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_S_uvw_fiPfj__, sz_);
#else
   d_S_uvw_fiPfj__ = clCreateBuffer(context__, CL_MEM_READ_ONLY, sz_, NULL, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_S_uvw_fiPfj__ device buffer");


   sz_ = sizeof(REAL);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_fi_abs__, sz_);
#else
   d_fi_abs__ = clCreateBuffer(context__, CL_MEM_READ_ONLY, sz_, NULL, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_fi_abs__ device buffer");
   
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_fj_abs__, sz_);
#else
   d_fj_abs__ = clCreateBuffer(context__, CL_MEM_READ_ONLY, sz_, NULL, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_fj_abs__ device buffer");
   
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_fiPfj_abs__, sz_);
#else
   d_fiPfj_abs__ = clCreateBuffer(context__, CL_MEM_READ_ONLY, sz_, NULL, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_fiPfj_abs__ device buffer");



#elif (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)


   sz_ = extdata__.NN__ * sizeof(REAL);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_nodal_vel__, sz_);
#else
   d_wind_nodal_vel__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_nodal_vel__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_nodal_vel__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_wind_nodal_vel__ on device");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_nodal_vel__, (void *)extdata__.wind_nodal_vel__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_wind_nodal_vel__.");
#endif


   sz_ = extdata__.NN__ * sizeof(int);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_nodal_windz__, sz_);
#else
   d_wind_nodal_windz__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_nodal_windz__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_nodal_windz__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_wind_nodal_windz__ on device");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_nodal_windz__, (void *)extdata__.wind_nodal_windz__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_wind_nodal_windz__.");
#endif


   sz_ = (size_t)(extdata__.NWZ__ * 3 * sizeof(REAL));
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_turb_std__, sz_);
#else
   d_wind_turb_std__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_turb_std__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_turb_std__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_wind_turb_std__ on device");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_turb_std__, (void *)extdata__.wind_turb_std__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_wind_turb_std__.");
#endif

   sz_ *= 3;
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_turb_scales__, sz_);
#else
   d_wind_turb_scales__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_turb_scales__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_turb_scales__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create d_wind_turb_scales__ on device");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_turb_scales__, (void *)extdata__.wind_turb_scales__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_wind_turb_scales__.");
#endif


# if (BSACL_KERNEL_ID==3)

   sz_ = nfi__ * sizeof(REAL);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_fi__, sz_);
#else
   d_fi__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)fi__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_fi__==NULL) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create   d_fi__   device buffer");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_fi__, (void *)fi__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_fi__.");
#endif


#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_fj__, sz_);
#else
   d_fj__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)fj__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create   d_fj__   device buffer");
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_fj__, (void *)fj__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to write device buffer  d_fj__.");
#endif


# endif // (BSACL_KERNEL_ID==3)


#endif // BSACL_KERNEL_ID



   if (extdata__.m3mf__ == NULL)
      ABORT_INTERNAL_RETURN_VOID(-83, "Result array was not acquired. Aborting.");
   sz_ = extdata__.DIM_M3_M__ * sizeof(REAL);

#if (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
   sz_ *= n_work_groups__[0];
#endif

#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_m3mf__, sz_);
#else
#  if (BSACL_KERNEL_ID==1)
   d_m3mf__ = clCreateBuffer(
      context__, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.m3mf__, &ierr_);
#  elif (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
   d_m3mf__ = clCreateBuffer(
      context__, CL_MEM_READ_WRITE, sz_, NULL, &ierr_);
#  endif
#endif
   if (ierr_ != BSACL_SUCCESS || d_m3mf__==NULL)
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to create  d_m3mf__  device buffer.");

#if (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
# ifndef BSACL_USE_CUDA__
   REAL zero_ = (REAL)0;
   ierr_ = clEnqueueFillBuffer(cqueues__[i_def_dev_id_], d_m3mf__, &zero_, 1, 0, sz_, 0, NULL, NULL);
   if (ierr_ != BSACL_SUCCESS)
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to enqueue fill of  d_m3mf__  device buffer.");

   ierr_ = clFinish(cqueues__[i_def_dev_id_]);
   if (ierr_ != BSACL_SUCCESS)
      ABORT_INTERNAL_RETURN_VOID(ierr_, 
         "Failed to finish default command queue when filling  d_m3mf__  device buffer.");
# else
   ierr_ = cudaMemset(d_m3mf__, 0, sz_);
   if (ierr_ != BSACL_SUCCESS)
      ABORT_INTERNAL_RETURN_VOID(ierr_, "Failed to enqueue fill of  d_m3mf__  device buffer.");
# endif // BSACL_USE_CUDA__

#endif

}








#ifndef BSACL_USE_CUDA__
size_t readFileContent_(char ** buf_, const char * fname_) {
   FILE *fp;
   errno_t nbytes_;
   long fsize_ = 0L;
   char ffp_[256];

   fsize_ = (long)(strlen(BASE_DIRECTORY) + strlen(fname_));
   strcpy(ffp_, BASE_DIRECTORY);
   nbytes_ = strcat_s(ffp_, fsize_ + 1, fname_);

   nbytes_ = fopen_s(&fp, ffp_, (const char *)"r");
   if (nbytes_ != 0) return 0;

   rewind(fp);
   fseek(fp, 0L, SEEK_END);
   fsize_ = ftell(fp);
   rewind(fp);

   *buf_ = (char *)malloc(sizeof(char) * (fsize_+1));
   if (*buf_ == NULL) return 0;

   nbytes_ = (int)fread_s(*buf_, fsize_, 1, fsize_, fp);
   if (nbytes_ == 0LL) {
      free(*buf_);
      ABORT_INTERNAL_RETURN_IERR(-431, "Failed to read source file.");
   }
   (*buf_)[nbytes_] = '\0'; // NULL terminate string

   return (size_t)nbytes_;
}




size_t assembleFinalCLSource_(char **clBaseStrings, char **evalFctStrings) {

   size_t fsize_ = 0, clsize_ = 0, evalsize_ = 0, fpos_ = 0;

   if (clBaseStrings != NULL)  clsize_   = strlen(*clBaseStrings);
   if (evalFctStrings != NULL) evalsize_ = strlen(*evalFctStrings);

   fsize_ = clsize_ + evalsize_;
   if (fsize_ == 0) return 0;

   fsize_ += 4;
   clSourceStrings__  = (char **)malloc(sizeof(char *));
   *clSourceStrings__ = (char  *)malloc(sizeof(char) * fsize_);
   strcpy(*clSourceStrings__, *clBaseStrings);
   fpos_ = clsize_;
   if (evalFctStrings && *evalFctStrings) {
      *(*clSourceStrings__ + fpos_) = '\n';
      fpos_++;
      *(*clSourceStrings__ + fpos_) = '\n';
      strcpy(*clSourceStrings__ + ++fpos_, *evalFctStrings);
      fpos_ += evalsize_+1;
   }
   *(*clSourceStrings__ + fpos_) = '\0';  // null terminate

// #ifdef BSACLC_DEBUG
//    printf("\n%s\n", *clSourceStrings__);
// #endif
   return fpos_;
}




void assembleProgramBuildOptsString_()
{
   char *buf = prog_compile_opts__;

   strcpy(buf, "-D BSACL_BASE_DIR=\"");
   buf += 19;
   strcpy(buf, BASE_DIRECTORY);
   buf += strlen(BASE_DIRECTORY);
   strcpy(buf, "\" \0");
   buf += 2;

   // BUG: user defined !
   strcpy(buf, "-D BSACL_CONV_PULSATION ");
   buf += 24;

   strcpy(buf, "-D BSACL_WIpWG=");
   buf += 15;
   strcpy(buf, STRINGIFYMACRO_VALUE(BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP));
   buf += strlen(STRINGIFYMACRO_VALUE(BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP));

   strcpy(buf, " -D BSACL_KERNEL_ID=");
   buf += 20;
   strcpy(buf, STRINGIFYMACRO_VALUE(BSACL_KERNEL_ID));
   ++buf;


#if (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
   strcpy(buf, " -D BSACL_WIND_PSD_ID=");
   buf += 22;
   strcpy(buf, STRINGIFYMACRO_VALUE(BSACL_PSD_TYPE_DAVENPORT));
   ++buf;

# ifdef BSACL_PASS_PARAMS_BY_MACRO
   strcpy(buf, " -D BSACL_PASS_PARAMS_BY_MACRO");
   buf += 30;

   strcpy(buf, " -D NTC__=");
   buf += 10;
   sprintf(buf, "%-1u", extdata__.NTC__);
   ++buf;

   strcpy(buf, " -D NNL__=");
   buf += 10;
   sprintf(buf, "%-10u", extdata__.NNODES_LOAD__);
   buf += 10;

   strcpy(buf, " -D NN__=");
   buf += 9;
   sprintf(buf, "%-10u", extdata__.NN__);
   buf += 10;

   strcpy(buf, " -D NM_EFF__=");
   buf += 13;
   sprintf(buf, "%-5u", extdata__.NMODES_EFF__);
   buf += 5;

// # if (BSACL_KERNEL_ID==3)
//    strcpy(buf, " -D NFI__=");
//    buf += 10;
//    sprintf(buf, "%-10u", nfi__);
//    buf += 10;

//    strcpy(buf, " -D NFJ__=");
//    buf += 10;
//    sprintf(buf, "%-10u", nfj__);
//    buf += 10;
// # endif

# endif // defined BSACL_PASS_PARAMS_BY_MACRO
#endif // (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)

   *buf = '\0';

#ifdef BSACLC_DEBUG
   printf("\n%s\n", prog_compile_opts__);
#endif
   return;
}





BSACL_INT buildCLProgram_() {
   BSACL_INT ierr_;
   size_t size_;
   char *clSource_;

   /** Create and build program. */
   size_ = readFileContent_(&clSource_, "bsacl.cl");
   size_ = assembleFinalCLSource_(&clSource_, evalFctStrings__);
   if (clSourceStrings__ == NULL) ABORT_INTERNAL_RETURN_IERR(-324, "Failed to assemble CL source strings.");
   free(clSource_), clSource_ = NULL;
   if (evalFctStrings__ != NULL) { free(*evalFctStrings__); evalFctStrings__ = NULL; }

   program__ = clCreateProgramWithSource(context__, 
      1, (const char **)clSourceStrings__, &size_, &ierr_);
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_IERR(ierr_, "Failed to create CL program object.");

   assembleProgramBuildOptsString_();
   ierr_ = clBuildProgram(program__, 1, devices__, prog_compile_opts__, NULL, NULL);
   if (ierr_ != BSACL_SUCCESS) {
      char buffer[20480];
      clGetProgramBuildInfo(program__, devices__[i_def_dev_id_],
         CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &size_);
      printf("\n Build log:\n\n%s\n", buffer);
      ABORT_INTERNAL_RETURN_IERR(ierr_, "Failed to build CL program.");
   }

#ifdef BSACL_GENERATE_PTX_BINARY
   FILE *ptxfile_ = NULL;
   char *buffer;
   
   size_  = 0llu;
   buffer = (char *)malloc(sizeof(char) * 1024 * 1000); // NOTE: need to allocate enough memory!
   ierr_  = clGetProgramInfo(program__, CL_PROGRAM_BINARIES, 20480, &buffer, &size_);
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_RETURN_IERR(ierr_, "Failed to query CL program binaries.");
   ptxfile_  = fopen("program.ptx", "w");
   if (ptxfile_ == NULL) ABORT_INTERNAL_RETURN_IERR(-544, "Failed to open file for writing PTX.");
   fprintf(ptxfile_, "%s", buffer);
   fclose(ptxfile_);
   free(buffer);
#endif

   return ierr_;
}





BSACL_INT setBfmKernelArgs_(void) 
{
   BSACL_INT ierr_  = 0;
   BSACL_UINT iarg_ = 0;

   REAL dInfl_ = (*(fi__ + 1) - *fi__) * 2 * BSACL_PI;  // rad/s
   dInfl_ = dInfl_ * dInfl_;  // infl area

#if (BSACL_KERNEL_ID==1) || (!defined(BSACL_PASS_PARAMS_BY_MACRO))
   ierr_  = clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(unsigned int), &extdata__.NTC__);
#endif
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem),       &d_tc__);
#if (BSACL_KERNEL_ID==1) || (!defined(BSACL_PASS_PARAMS_BY_MACRO))
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(unsigned int), &extdata__.NNODES_LOAD__);
#endif
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem),       &d_nodes_load__);
   
#if (BSACL_KERNEL_ID==1)
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem), &d_S_uvw_fi__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem), &d_fi_abs__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem), &d_S_uvw_fj__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem), &d_fj_abs__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem), &d_S_uvw_fiPfj__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem), &d_fiPfj_abs__);
#elif (BSACL_KERNEL_ID==2)
   iarg_ = 6;
#elif (BSACL_KERNEL_ID==3)
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem),       &d_fi__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(unsigned int), &nfi__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(cl_mem),       &d_fj__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++,  sizeof(unsigned int), &nfj__);
#endif
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(REAL),   &dInfl_);

#if (BSACL_KERNEL_ID==1) || (!defined(BSACL_PASS_PARAMS_BY_MACRO))
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(unsigned int), &extdata__.NMODES_EFF__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(unsigned int), &extdata__.NDEGW__);
#endif
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(cl_mem),       &d_phiTc__);
#if (BSACL_KERNEL_ID==1) || (!defined(BSACL_PASS_PARAMS_BY_MACRO))
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(unsigned int), &extdata__.NN__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(unsigned int), &extdata__.NNOD_CORR__);
#endif
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(cl_mem),       &d_nod_corr__);

#if (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(cl_mem), &d_wind_nodal_vel__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(cl_mem), &d_wind_turb_scales__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(cl_mem), &d_wind_turb_std__);
   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(cl_mem), &d_wind_nodal_windz__);
#endif

   ierr_ |= clSetKernelArg(kernel_bfm__, iarg_++, sizeof(cl_mem), &d_m3mf__);
   return ierr_;
}
#endif // BSACL_USE_CUDA__







ierr_t getOptKernelDims_() {
   size_t ntwi_, ntot_;
   BSACL_UINT nwg_;

   // BUG: use power function ??
   ntwi_ = extdata__.NMODES_EFF__  * extdata__.NMODES_EFF__  * extdata__.NMODES_EFF__;
   ntot_ = extdata__.NNODES_LOAD__ * extdata__.NNODES_LOAD__ * extdata__.NNODES_LOAD__ * ntwi_;

#if (BSACL_KERNEL_ID==1)
   // get theoretical n. of WG using max WIpWG
# ifdef BSACL_USE_CUDA__
   nwg_ = BSACL_OPT_N_WORK_GROUPS;
# else
   nwg_ = (BSACL_UINT)ceil((double)ntwi_ / dev_max_work_group_size__[i_def_dev_id_]);
   if (nwg_ < BSACL_OPT_N_WORK_GROUPS)
      nwg_ = (BSACL_UINT)BSACL_OPT_N_WORK_GROUPS;
# endif

   size_t wi_p_wg_ = (size_t)ceil((double)ntwi_ / nwg_);
   if (wi_p_wg_ < BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP) {
      wi_p_wg_ = BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP;
      nwg_     = (BSACL_UINT)ceil((double)ntwi_ / wi_p_wg_);
   }

   n_work_dims__ = 1;
   n_work_groups__[0] = nwg_;
   n_work_groups__[1] = 1;
   n_work_groups__[2] = 1;

   global_dims_ie_NTWI__[0] = wi_p_wg_ * nwg_;
   global_dims_ie_NTWI__[1] = 1;
   global_dims_ie_NTWI__[2] = 1;

   local_dims_ie_WIpWG__[0] = wi_p_wg_;
   local_dims_ie_WIpWG__[1] = 1;
   local_dims_ie_WIpWG__[2] = 1;

#elif (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)

   n_work_dims__ = 2;

   global_dims_ie_NTWI__[1] = ntwi_; // NM^3

   ntwi_ = extdata__.NNODES_LOAD__ * extdata__.NNODES_LOAD__ * extdata__.NNODES_LOAD__;
   nwg_  = (BSACL_UINT)ceil((double)ntwi_ / BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP);

   n_work_groups__[0] = nwg_;
   n_work_groups__[1] = global_dims_ie_NTWI__[1];
   
   global_dims_ie_NTWI__[0] = BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP * nwg_;

   local_dims_ie_WIpWG__[0] = BSACL_OPT_N_WORK_ITEMS_PER_WORK_GROUP;
   local_dims_ie_WIpWG__[1] = 1;

# ifdef BSACL_USE_CUDA__
   cuda_nblocks__.x           = (unsigned int)n_work_groups__[0];
   cuda_nblocks__.y           = (unsigned int)n_work_groups__[1];
   cuda_nblocks__.z           = 1;
   cuda_nthreads_perblock__.x = (unsigned int)local_dims_ie_WIpWG__[0];
   cuda_nthreads_perblock__.y = (unsigned int)local_dims_ie_WIpWG__[1];
   cuda_nthreads_perblock__.z = 1;
# endif

#endif 


   ntot_ *= extdata__.NTC__;


#ifdef BSACLC_DEBUG
   printf("\n\n%sTotal n. of operations (%u^3 x %u^3 x %u)  %llu  (per freqs. pair).", 
      INFO_MSG, extdata__.NNODES_LOAD__, extdata__.NMODES_EFF__, extdata__.NTC__, ntot_);
   printf("\n%sEnqueueing kernel using:", INFO_MSG);
   printf("\n%s - n. work dimensions : %u", CONT_MSG, n_work_dims__);
   printf("\n%s - n. work-groups     : %u", CONT_MSG, nwg_);
   printf("\n%s - n. work-items / WG : %llu - %llu - %llu", 
      CONT_MSG, local_dims_ie_WIpWG__[0], local_dims_ie_WIpWG__[1], local_dims_ie_WIpWG__[2]);
   printf("\n%s - n. tot work-items  : %llu - %llu - %llu\n\n", 
      CONT_MSG, global_dims_ie_NTWI__[0], global_dims_ie_NTWI__[1], global_dims_ie_NTWI__[2]);
#endif

#if (BSACL_KERNEL_ID==1)
   if (nwg_ < MIN_N_OF_WORK_GROUPS) return (ierr_t)BSACL_PROBLEM_DIMENSIONS_TOO_SMALL_;
#endif

   return BSACL_SUCCESS;
}








//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//
//                        SECTION:   EXPOSED library functions (bsacl.h)
//
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------



/**
 * @brief Initialises CL runtime
 * @note  External data has to be fully initialised for this function to work!
 * 
 * @return (cl_int) Error code
 */
void bsaclInit(int *__EXT_PTR_CONST ierr)
{

   if (has_halted__ == 1U) { *ierr = -999; return; }

#ifdef BSACL_USE_CUDA__
   *ierr = BSACL_SUCCESS;
   return;
#else
   BSACL_INT ierr_;

   // Acquire platform
   ierr_ = getCLPlatformFromInfo_(&platform__, CL_PLATFORM_NAME, "NVIDIA\0", 0);
   if (ierr_ != BSACL_SUCCESS) goto ret_;

   // Query devices by type
   n_devices__ = getCLDevicesByType_(&platform__, device_type__, &devices__);
   if (n_devices__ <= 0) {
      ierr_ = (cl_int) n_devices__;
      goto ret_;
   }

   getCLDevicesInfo_();

   // get CL context
   context__ = clCreateContext(NULL, n_devices__, devices__, NULL, NULL, &ierr_);
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create context", ret_);

   // One command queue per device (for the moment)
   cqueues__ = (cl_command_queue *) malloc(n_devices__ * sizeof(cl_command_queue));
   if (cqueues__ == NULL) {
      ierr_ = (cl_int)BSACL_CQUEUES_CREATION_ERROR_;
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create command queues.", ret_);
   }
   for (unsigned i = 0; i < n_devices__; ++i) {
      cqueues__[i] = clCreateCommandQueue(context__, devices__[i], NULL, &ierr_);
      if (ierr_ != BSACL_SUCCESS) {
         printf("\n\n%sFailed to create  %u-th  command queue. Aborting (%d).", ERRR_MSG, i, ierr_);
         abortInternal_(0, NULL);
         goto ret_;
      }
   }

   // NOTE: put it here so that we avoid set to 1 (true) if error occurs
   bsacl_init_called__ = 1U;

   ret_: *ierr = (int)ierr_;
   return;
#endif
}




/**
 * @brief BSACL core function. Enqueue selected kernel to run on GPU(s).
 * 
 * @param ierr Error code.
 */
void bsaclRun(int *__EXT_PTR_CONST ierr) {
   ierr_t ierr_;
#ifndef BSACL_USE_CUDA__
   char * clSource_ = NULL;
#endif

   const size_t n_freqs_tot_ = nfj__ * nfi__;
   size_t size_ = 0ull;

#if (BSACL_KERNEL_ID!=3)
   REAL fi_, fj_;
   size_t icount_ = 0ull;
#endif
#if (BSACL_KERNEL_ID==1)
   REAL fi_abs_, fj_abs_, fiPfj_, fiPfj_abs_;
#endif

#ifndef BSACL_USE_CUDA__
   cl_event mainKerEv_;
#elif (defined BSACL_USE_CUDA__) && (BSACL_KERNEL_ID!=1)
   REAL dInfl_ = (*(fi__ + 1) - *fi__) * 2 * BSACL_PI;  // rad/s
   dInfl_ = dInfl_ * dInfl_;  // infl area
#endif


   if (has_halted__ == 1U) {
      ierr_ = (ierr_t)-1;
      goto ret_;
   }

   if (fi__ == NULL || fj__ == NULL) 
      { ierr_ = (ierr_t)-20; printMsg_(ERRR_MSG, "No computation frequencies acquired."); goto ret_; }

   ierr_ = BSACL_SUCCESS;
   if (bsacl_init_called__ == 0U) bsaclInit((int *)&ierr_);
   if (ierr_ != BSACL_SUCCESS) { printMsg_(ERRR_MSG, "Failed to initialise BSACL."); goto ret_; }

   /** Define Kernel dimensions */
   ierr_ = getOptKernelDims_();
   if (ierr_ != BSACL_SUCCESS) goto ret_;


#ifndef BSACL_USE_CUDA__
   ierr_ = buildCLProgram_();
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create CL program object.", ret_);

   /** Create kernel */
   kernel_bfm__ = clCreateKernel(program__, "bfm_kernel", &ierr_);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create CL kernel object.", ret_);
#endif // not def BSACL_USE_CUDA__
   

#if (BSACL_KERNEL_ID==1)
   if (BASE_PSD_ALLOC_SIZE__ == 0U)
      BASE_PSD_ALLOC_SIZE__ = extdata__.NNODES_LOAD__ * extdata__.NTC__;
   size_ = BASE_PSD_ALLOC_SIZE__ * sizeof(REAL);
#endif


   initCreateDeviceBuffers_();


#ifndef BSACL_USE_CUDA__
   ierr_ = setBfmKernelArgs_();
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to set kernel arguments (1).", ret_);
#endif // not def BSACL_USE_CUDA__



   // ----------------------------
   //            MAIN LOOP
   // ----------------------------

   printf("\n\n%sReady to perform  %12llu  iterations.\n", INFO_MSG, n_freqs_tot_);

#if (BSACL_KERNEL_ID!=3)
   // for (unsigned int j_ = 0; j_ < 1; ++j_) {
   for (unsigned int j_ = 0; j_ < nfj__; ++j_) {

      fj_ = *(fj__ + j_);
#endif

#if (BSACL_KERNEL_ID==1)

#ifdef BSACL_USE_CUDA__
      ierr_ = cudaMemcpy(d_S_uvw_fj__, (S_uvw__ + (j_ * BASE_PSD_ALLOC_SIZE__)), size_, cudaMemcpyHostToDevice);
#else
      ierr_ = clEnqueueWriteBuffer(
            cqueues__[i_def_dev_id_], d_S_uvw_fj__, CL_FALSE, 0, size_, 
               (S_uvw__ + (j_ * BASE_PSD_ALLOC_SIZE__)), 0, NULL, NULL);
#endif
      if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_S_uvw_fj__.", ret_);


      fj_abs_ = fabs(fj_);
#ifdef BSACL_USE_CUDA__
      ierr_ = cudaMemcpy(d_fj_abs__, &fj_abs_, sizeof(REAL), cudaMemcpyHostToDevice);
#else
      ierr_   = clEnqueueWriteBuffer(
            cqueues__[i_def_dev_id_], d_fj_abs__, CL_FALSE, 0, sizeof(REAL), 
               &fj_abs_, 0, NULL, NULL);
#endif
      if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_fj_abs__.", ret_);

#elif (BSACL_KERNEL_ID==2)

      ierr_ = clSetKernelArg(kernel_bfm__, 5U, sizeof(REAL), &fj_);
      if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to set   fj_   kernel argument.", ret_);

#endif // (BSACL_KERNEL_ID)


#if (BSACL_KERNEL_ID!=3)
      // for (unsigned int i_ = 0; i_ < 1; ++i_) {
      for (unsigned int i_ = 0; i_ < nfi__; ++i_) {

         fi_ = *(fi__ + i_);
#endif

#if (BSACL_KERNEL_ID==1)

#ifdef BSACL_USE_CUDA__
         ierr_ = cudaMemcpy(d_S_uvw_fi__, (S_uvw__ + (i_ * BASE_PSD_ALLOC_SIZE__)), size_, cudaMemcpyHostToDevice);
#else
         ierr_ = clEnqueueWriteBuffer(
            cqueues__[i_def_dev_id_], d_S_uvw_fi__, CL_FALSE, 0, size_, 
               (S_uvw__ + (i_ * BASE_PSD_ALLOC_SIZE__)), 0, NULL, NULL);
#endif
         if (ierr_ != BSACL_SUCCESS) 
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_S_uvw_fi__.", ret_);

         
         fi_abs_ = fabs(fi_);
#ifdef BSACL_USE_CUDA__
         ierr_ = cudaMemcpy(d_fi_abs__, &fi_abs_, sizeof(REAL), cudaMemcpyHostToDevice);
#else
         ierr_   = clEnqueueWriteBuffer(
               cqueues__[i_def_dev_id_], d_fi_abs__, CL_FALSE, 0, sizeof(REAL), 
                  &fi_abs_, 0, NULL, NULL);
#endif
         if (ierr_ != BSACL_SUCCESS) 
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_fi_abs__.", ret_);

         // Evaluate PSD at (fi+fj)
         fiPfj_ = fi_ + fj_;
         for (unsigned itc = 0; itc < extdata__.NTC__; ++itc)
            evalFunc_ptr__(
               extdata__.tc__[itc], 1, &fiPfj_, extdata__.NNODES_LOAD__, 
                  (S_uvw_fiPfj__ + (itc * extdata__.NNODES_LOAD__)));

#ifdef BSACL_USE_CUDA__
         ierr_ = cudaMemcpy(d_S_uvw_fiPfj__, S_uvw_fiPfj__, size_, cudaMemcpyHostToDevice);
#else
         ierr_ = clEnqueueWriteBuffer(
            cqueues__[i_def_dev_id_], d_S_uvw_fiPfj__, CL_FALSE, 0, size_, 
               S_uvw_fiPfj__, 0, NULL, NULL);
#endif
         if (ierr_ != BSACL_SUCCESS) 
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_S_uvw_fiPfj__.", ret_);


         fiPfj_abs_ = fabs(fiPfj_);
#ifdef BSACL_USE_CUDA__
         ierr_ = cudaMemcpy(d_fiPfj_abs__, &fiPfj_abs_, sizeof(REAL), cudaMemcpyHostToDevice);
#else
         ierr_ = clEnqueueWriteBuffer(
               cqueues__[i_def_dev_id_], d_fiPfj_abs__, CL_FALSE, 0, sizeof(REAL), 
                  &fiPfj_abs_, 0, NULL, NULL);
#endif
         if (ierr_ != BSACL_SUCCESS) 
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_fiPfj_abs__.", ret_);

#elif (BSACL_KERNEL_ID==2)

         ierr_ = clSetKernelArg(kernel_bfm__, 4U, sizeof(REAL), &fi_);
         if (ierr_ != BSACL_SUCCESS) 
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to set   fi_   kernel argument.", ret_);
#endif // (BSACL_KERNEL_ID==1)


#ifdef BSACL_USE_CUDA__
         bfm_kernel<<<cuda_nblocks__, cuda_nthreads_perblock__>>>(\
            extdata__.PSD_ID__,
            extdata__.NTC__, 
            d_tc__,
            extdata__.NNODES_LOAD__,
            d_nodes_load__,
#if (BSACL_KERNEL_ID==1)
            d_S_uvw_fi__,
            d_fi_abs__,
            d_S_uvw_fj__,
            d_fj_abs__,
            d_S_uvw_fiPfj__,
            d_fiPfj_abs__,
#elif (BSACL_KERNEL_ID==2)
            fi_,
            fj_,
#elif (BSACL_KERNEL_ID==3)
            d_fi__,
            nfi__,
            d_fj__,
            nfj__,
#endif
            dInfl_,
            extdata__.NMODES_EFF__,
            extdata__.NDEGW__,
            d_phiTc__,
            extdata__.NN__,
            extdata__.NNOD_CORR__,
            d_nod_corr__,
#if (BSACL_KERNEL_ID!=1)
            d_wind_nodal_vel__,
            d_wind_turb_scales__,
            d_wind_turb_std__,
            d_wind_nodal_windz__,
#endif
            d_m3mf__);
         ierr_ = cudaDeviceSynchronize();
         // ierr_ = cudaGetLastError();
#else
         ierr_ = clEnqueueNDRangeKernel(
            cqueues__[i_def_dev_id_], kernel_bfm__, n_work_dims__, NULL, 
               global_dims_ie_NTWI__, local_dims_ie_WIpWG__, 
                  0, NULL, &mainKerEv_);
#endif // BSACL_USE_CUDA__
         if (ierr_ != BSACL_SUCCESS) 
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to enqueue kernel execution.", ret_);


#ifndef BSACL_USE_CUDA__
         clWaitForEvents(1, &mainKerEv_);
         // ierr_ = clFinish(cqueues__[i_def_dev_id_]);
         // if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to finish command queue (2).", ret_);
#endif


#if (BSACL_KERNEL_ID!=3)
         ++icount_;
      
      } // nfi__


      printf("\n%s  %12llu  out of  %12llu  iterations done...", 
         INFO_MSG, icount_, n_freqs_tot_);
   } // nfj__
#endif

   printf("\n\n");

   // clReleaseEvent(mainKerEv_);
#if (BSACL_KERNEL_ID==1)

# ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy(extdata__.m3mf__, d_m3mf__, extdata__.DIM_M3_M__ * sizeof(REAL), cudaMemcpyDeviceToHost);
# else
   ierr_  = clEnqueueReadBuffer(cqueues__[i_def_dev_id_], d_m3mf__, CL_TRUE, 
         0, extdata__.DIM_M3_M__ * sizeof(REAL), extdata__.m3mf__, 0, NULL, NULL);
# endif

#elif (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
   
   REAL *rtmp_;
   size_ = sizeof(REAL) * extdata__.DIM_M3_M__*n_work_groups__[0];
   rtmp_ = (REAL *)malloc(size_);
# ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy(rtmp_, d_m3mf__, size_, cudaMemcpyDeviceToHost);
# else
   ierr_  = clEnqueueReadBuffer(cqueues__[i_def_dev_id_], d_m3mf__, CL_TRUE, 0, size_, rtmp_, 0, NULL, NULL);
# endif

#endif // BSACL_KERNEL_ID
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to read results from device.", ret_);


#if (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)
   /** Cumulate all single WG contributions */
   for (BSACL_UINT iwgx_ = 0; iwgx_ < n_work_groups__[0]; iwgx_++) {
      for (BSACL_UINT i_ = 0; i_ < extdata__.DIM_M3_M__; i_++) {
         extdata__.m3mf__[i_] += rtmp_[iwgx_*extdata__.DIM_M3_M__ + i_];
      }
   }
   free(rtmp_);
   rtmp_ = NULL;
#endif


   ierr_ = releaseDeviceMem_();
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to release device memory.", ret_);

#ifndef BSACL_USE_CUDA__
   ierr_ = clFinish(cqueues__[i_def_dev_id_]);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to finish default command queue.", ret_);

   ierr_ = clReleaseCommandQueue(cqueues__[i_def_dev_id_]);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to release command queue.", ret_);

   ierr_ = clReleaseKernel(kernel_bfm__);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to release kernel.", ret_);

   ierr_ = clReleaseProgram(program__);
   if (ierr_ != BSACL_SUCCESS) 
      abortInternal_(ierr_, "Failed to release program.");
#endif

ret_: *ierr = (int)ierr_;
#ifndef BSACL_USE_CUDA__
   if (clSource_) free(clSource_);
#endif
   return;
}











void bsaclAcquirePSDId(const uint32_t psdid) {
   if (has_halted__ == 1U) return;
   extdata__.PSD_ID__ = psdid;
}




/**
 * @brief Get pointer to memory containing Modal Matrix
 * 
 * @note For dimensions, no need to pass any reference. Copy.
 * 
 * @param modmat pointer to base memory address
 * @param natf   vector of modal natural frequencies
 * @param ndofs  n. of total degrees of freedom
 * @param nmodes n. of modes (used)
 */
void bsaclAcquireStructModMat(
      REAL *__EXT_PTR_CONST modmat, REAL *__EXT_PTR_CONST natf, const uint32_t ndofs, const uint32_t nmodes) {

   if (has_halted__ == 1U) return;
   extdata__.modmat__   = modmat;
   extdata__.natfreqs__ = natf;
   extdata__.NDOFS__    = ndofs;
   extdata__.NMODES__   = nmodes;
// #ifdef BSACLC_DEBUG
//    for (unsigned i = 0; i < ndofs; ++i) {
//       printf("\n");
//       for (unsigned j = 0; j < nmodes; ++j) {
//          printf("%10.4g  ", extdata__.modmat__[j*ndofs + i]);
//       }
//    } 
// #endif
}


/**
 * @brief Get ptr to memory containing loaded nodes list
 * 
 * @param nodes_load 
 * @param nnodes_l 
 */
void bsaclAcquireLoadedNodesList(int *__EXT_PTR_CONST nodes_load, const uint32_t nnodes_l) {

   if (has_halted__ == 1U) return;
   if (extdata__.NNODES_LOAD__ != 0 && nnodes_l != extdata__.NNODES_LOAD__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of loaded nodes does not match.");
   extdata__.nodes_load__  = nodes_load;
   extdata__.NNODES_LOAD__ = nnodes_l;
// #ifdef BSACLC_DEBUG
//    printf("\n");
//    for (unsigned i = 0; i < nnodes_l; ++i) printf("  %4d", extdata__.nodes_load__[i]);
// #endif
}


void bsaclAcquireTotalNOfNodes(const uint32_t nnodes_tot) {
   if (has_halted__ == 1U) return;
   if (extdata__.NN__ != 0 && nnodes_tot != extdata__.NN__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of total nodes does not match.");
   extdata__.NN__  = nnodes_tot;
}


/**
 * @brief Get ptr to memory holding list of effective used modes
 * 
 * @param modes 
 * @param nmodes_eff 
 */
void bsaclAcquireUsedModesList(int *__EXT_PTR_CONST modes, const uint32_t nmodes_eff) {

   if (has_halted__ == 1U) return;
   if (extdata__.NMODES_EFF__ != 0 && nmodes_eff != extdata__.NMODES_EFF__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of effective modes does not match.");
   extdata__.modes_eff__  = modes;
   extdata__.NMODES_EFF__ = nmodes_eff;
// #ifdef BSACLC_DEBUG
//    printf("\n");
//    for (unsigned i = 0; i < nmodes_eff; ++i) printf("  %4d", extdata__.modes_eff__[i]);
// #endif
}


/**
 * @brief Get ptr to memory holding Wind Force Coefficients data
 * 
 * @param wfc       wind force coefficient matrix
 * @param nnodes_l  n. of loaded nodes
 * @param nlibs     n. of degrees of freedom per node
 * @param ndegw     n. of coefficients per element (transformation degree)
 */
void bsaclAcquireWindCoeffs(REAL *__EXT_PTR_CONST wfc, const uint32_t nnodes_l, const uint32_t nlibs, const uint32_t ndegw) {

   if (has_halted__ == 1U) return;
   if (extdata__.NNODES_LOAD__ != 0 && nnodes_l != extdata__.NNODES_LOAD__) 
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of loaded nodes does not match.");
   if (extdata__.NLIBS__ != 0 && nlibs != extdata__.NLIBS__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of LIBS does not match.");
   if (extdata__.NDEGW__ != 0 && ndegw != extdata__.NDEGW__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "DEGW does not match.");
   extdata__.NNODES_LOAD__ = nnodes_l;
   extdata__.NLIBS__       = nlibs;
   extdata__.NDEGW__       = ndegw;
   extdata__.wfc__         = wfc;
}


/**
 * @brief Get ptr to memory holding list of effective turbulent components considered
 *
 * @param tc   list of turbulent components
 * @param ntc  n. of turb. components in the list
 */
void bsaclAcquireTurbComponentsList(int *__EXT_PTR_CONST tc, const uint32_t ntc) {

   if (has_halted__ == 1U) return;
   extdata__.NTC__ = ntc;
   extdata__.tc__  = tc;
#if (BSACL_KERNEL_ID==1)
   if (extdata__.NNODES_LOAD__ != 0) BASE_PSD_ALLOC_SIZE__ = extdata__.NNODES_LOAD__ * ntc;
#endif
}


/**
 * @brief Get ptr to memory holding the product of Phi x C
 * 
 * @param phi_T_c    product of modal matrix times wind-force-coeffs matrix
 * @param nnodes_l   n. of loaded nodes
 * @param nmodes_eff n. of effective modes
 * @param ndegw      n. of coefficients per element (transformation degree)
 */
void bsaclAcquirePhiTimesCMat(REAL *__EXT_PTR_CONST phi_T_c, const uint32_t nmodes_eff, const uint32_t nnodes_l, const uint32_t ndegw) {

   if (has_halted__ == 1U) return;
   if (extdata__.NNODES_LOAD__ != 0 && nnodes_l != extdata__.NNODES_LOAD__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of loaded nodes does not match.");
   if (extdata__.NMODES_EFF__ != 0 && nmodes_eff != extdata__.NMODES_EFF__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of effective modes does not match.");
   if (extdata__.NDEGW__ != 0 && ndegw != extdata__.NDEGW__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "DEGW does not match.");
   extdata__.NNODES_LOAD__ = nnodes_l;
   extdata__.NMODES_EFF__  = nmodes_eff;
   extdata__.NDEGW__       = ndegw;
   extdata__.phi_T_c__     = phi_T_c;
}




void bsaclAcquireNodalCorrelation(REAL *__EXT_PTR_CONST nod_corr, const uint32_t nnod_corr) {
   if (has_halted__ == 1U) return;
   extdata__.nod_corr__  = nod_corr;
   extdata__.NNOD_CORR__ = nnod_corr;
}



void bsaclAcquireWindNodalVelocities(REAL *__EXT_PTR_CONST nod_vel) {
   if (has_halted__ == 1U) return;
   extdata__.wind_nodal_vel__ = nod_vel;
}

void bsaclAcquireWindNodalWindZones(int *__EXT_PTR_CONST nod_wz) {
   if (has_halted__ == 1U) return;
   extdata__.wind_nodal_windz__ = nod_wz;
}

void bsaclAcquireWindTurbScales(REAL *__EXT_PTR_CONST wt_scl, const uint32_t nwz) {
   if (has_halted__ == 1U) return;
   extdata__.NWZ__ = nwz;
   extdata__.wind_turb_scales__ = wt_scl;
}

void bsaclAcquireWindTurbStd(REAL *__EXT_PTR_CONST wt_std, const uint32_t nwz) {
   if (has_halted__ == 1U) return;
   extdata__.NWZ__ = nwz;
   extdata__.wind_turb_std__ = wt_std;
}





void bsaclAcquireEvalFunc(evalFct_t fct) {

   if (has_halted__ == 1U) return;
#if (BSACL_KERNEL_ID==1)
   if (fct == NULL) ABORT_INTERNAL_RETURN_VOID(-11, "Passed evaluation function pointer is NULL.");
   evalFunc_ptr__ = fct;
#endif
   return;
}


void bsaclAcquireEvalFuncByStrings(char **__EXT_PTR_CONST strings) {
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   size_t ilen_ = strlen(*strings);
   evalFctStrings__  = (char **)malloc(sizeof(char *));
   *evalFctStrings__ = (char  *)malloc(sizeof(char) * (ilen_ + 1));
   strcpy(*evalFctStrings__, *strings);
   *(*evalFctStrings__ + ilen_) = '\0';
#endif
   return;
}


void bsaclAcquireEvalFuncByFile(char *__EXT_PTR_CONST filename) {
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   size_t ierr_   = 0;
   char *fstrings = NULL;
   ierr_ = readFileContent_(&fstrings, filename);
   bsaclAcquireEvalFuncByStrings(&fstrings);
   free(fstrings);
#endif
   return;
}





void bsaclAcquireComputationFreqs(const uint32_t nfi, REAL *__EXT_PTR_CONST fi, const uint32_t nfj, REAL *__EXT_PTR_CONST fj) {
   if (has_halted__ == 1U) return;
   nfi__ = nfi;
   fi__  = fi;
   nfj__ = nfj;
   fj__  = fj;
   return;
}


void bsaclAcquireBaseWindTurbPSD(REAL *__EXT_PTR_CONST S_uvw) {
   if (has_halted__ == 1U) return;
   S_uvw__ = S_uvw;
   return;
}


void bsaclAcquireResultBFMVect(REAL *__EXT_PTR_CONST m3mf, const uint32_t idim) {
   if (has_halted__ == 1U) return;
   extdata__.m3mf__     = m3mf;
   extdata__.DIM_M3_M__ = idim;
   return;
}




void bsaclSetDeviceType(const uint32_t itype) {
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   switch (itype) {
      case BSACL_DEVICE_TYPE_CPU_:
         device_type__ = CL_DEVICE_TYPE_CPU;
         break;
      case BSACL_DEVICE_TYPE_GPU_:
         device_type__ = CL_DEVICE_TYPE_GPU;
         break;
      case BSACL_DEVICE_TYPE_ACC_:
         device_type__ = CL_DEVICE_TYPE_ACCELERATOR;
         break;
      case BSACL_DEVICE_TYPE_DEF_:
         break;
      default:
         abortInternal_(BSACL_INVALID_DEVICE_TYPE_, "Invalid device type requested.");
   }
#endif
   return;
}





void bsaclVerifyMaxAllocCondition(size_t idim, unsigned int *__EXT_PTR_CONST ican) {
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   if (dev_max_mem_alloc_size__[0] < idim * sizeof(REAL)) {
      *ican = 0U;
      return;
   }
#endif
   *ican = 1U;
   return;
}




void bsaclAbort(const int ierr) {
   abortInternal_(ierr, NULL);
}

/**
 * @brief Clens BSACL Runtime memory 
 * 
 */
void bsaclFinalise(void) {
   freeMem_();
}






// extern "C" {
#ifdef __cplusplus
}
#endif

