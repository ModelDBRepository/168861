//
//
// File author(s):  <Gabriel Rucabado and Antonio Garcia Dopico>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef _CUDAEXCEPTION
#define _CUDAEXCEPTION

#include <string>
#include <sstream>
#include <exception>
// Includes from cuda
#include "cuda_runtime.h"
#include <cublas_v2.h>
#include <cusparse_v2.h>

//
// DEFINITIONS
//
#define CUDA_ID_DEF 0

#define OK_CUSPARSE CUSPARSE_STATUS_SUCCESS
#define OK_CUBLAS CUBLAS_STATUS_SUCCESS
#define OK_CUDA cudaSuccess


//#######################
//#### CONFIGURATION ####
//#######################
#define NUM_TH_XAXPY 512      //<=1024 (&& ^2 with VER32)
#define NUM_TH_REDUCE 256     //<=1024 (&& ^2 with VER32)
#define NUM_TH_VECTOR_OP 512  //<=1024 (&& ^2 with VER32)
#define NUM_TH_MATRIX_OP 512  //<=1024 (&& ^2 with VER32)
#define NUM_LEC_TH 4          // 
#define MIN_MAJOR 2           // xa multiGPU
#define MIN_MINOR 0           // xa multiGPU

#define TH_HH 25
#define TH_BRANCH 32
#define TH_TERM 2
#define TH_ALL 64



//
// MACROS
//
#define DIM3(var, Vx, Vy, Vz)			\
  var.x = Vx;					\
  var.y = Vy;					\
  var.z = Vz						

#if DEBUG >= DEBUG_MESS_LVL
# define _DEBUG_CUDA(fmt, args...)					\
  printf (DEBUG_COLOR "DEBUG_CUDA: Grid<%05d-%05d-%05d> Block<%05d-%05d-%05d>" \
	  " [%s -.- %s:%05d] ::: "fmt DEFAULT"\n",			\
	  blockIdx.x, blockIdx.y, blockIdx.z,				\
	  threadIdx.x, threadIdx.y, threadIdx.z,			\
	  __BASE_FILE__, __FUNCTION__, __LINE__, ##args);

#else
# define _DEBUG_CUDA(d,f...)
#endif

# define _INFO_CUDA(fmt, args...)					\
  printf (INFO_COLOR "INFO_CUDA: Grid<%05d-%05d-%05d> Block<%05d-%05d-%05d>" \
	  " [%s -.- %s:%05d] ::: "fmt DEFAULT "\n",			\
	  blockIdx.x, blockIdx.y, blockIdx.z,				\
	  threadIdx.x, threadIdx.y, threadIdx.z,			\
	  __BASE_FILE__, __FUNCTION__, __LINE__, ##args);
# define _ERROR_CUDA(fmt, args...)					\
  fprintf (stderr, ERROR_COLOR "ERROR_CUDA: Grid<%05d-%05d-%05d> Block<%05d-%05d-%05d>" \
	   " [%s -.- %s:%05d] ::: "fmt DEFAULT"\n",			\
	   blockIdx.x, blockIdx.y, blockIdx.z,				\
	   threadIdx.x, threadIdx.y, threadIdx.z,			\
	   __BASE_FILE__, __FUNCTION__, __LINE__, ##args); 

# define CUDA_ALLOC(ptr, size)						\
  {									\
    cudaError_t error;							\
    std::stringstream params;						\
    if (( error = cudaMalloc ((void**)ptr, size)) != OK_CUDA){		\
      params << "ptr=" << ptr << " size=" << size;			\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,	\
			  "When allocating the memory", error, params.str());	\
    }									\
  }

#if DEBUG >= DEBUG_KERNEL_LVL
#define CUDA_EVENT_RECORD(ev, st)					\
  {									\
    cudaError_t error;							\
    if ((error = cudaEventRecord(ev, st)) != OK_CUDA){			\
      std::stringstream params;						\
      params << "event=" << ev << " stream=" << st;			\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,	\
			  "Tracking event updateHH<->updateMat",	\
			  error, params.str());				\
    }									\
  }

#define CUDA_EVENT_WAIT(ev, st)						\
  {									\
    cudaError_t error;							\
    if ((error = cudaStreamWaitEvent (st, ev, 0)) != OK_CUDA){		\
      std::stringstream params;						\
      params << "event=" << ev << " stream=" << st;			\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,	\
			  "Synchronizing updateHH<->updateMat",		\
			  error,  params.str());			\
    }									\
  }

#define CUDA_EVENT_SYNC(ev)						\
  {									\
    cudaError_t error;							\
    if ((error = cudaEventSynchronize (ev)) != OK_CUDA){		\
      std::stringstream params;						\
      params << "event=" << ev;					\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,	\
			  "Synchronizing updateHH<->updateMat",		\
			  error,  params.str());			\
    }									\
  }
#else
#define CUDA_EVENT_RECORD(ev, st)					\
  {									\
   cudaEventRecord(ev, st);						\
  }
#define CUDA_EVENT_WAIT(ev, st)						\
  {									\
   cudaStreamWaitEvent (st, ev, 0);					\
  }
#define CUDA_EVENT_SYNC(ev)						\
  {									\
    cudaEventSynchronize (ev);						\
  }
#endif


#define CUDA_CPY_SYNC(tar, src, len)					\
  {									\
    cudaError_t error;							\
    std::stringstream params;						\
    error = cudaMemcpy (tar, src, len, cudaMemcpyDefault);		\
    if (error != OK_CUDA){						\
      params << "src=" << src << " tar=" << tar << " len=" << len;	\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,	\
			  "Doing the memory transfer",	\
			  error, params.str());				\
    }									\
  }

#define CUDA_CPY(tar, src, len, st)					\
  {									\
    cudaError_t error;							\
    std::stringstream params;						\
    if ((error = cudaMemcpyAsync(tar, src, len, cudaMemcpyDefault, st)) != OK_CUDA){ \
      params << "src=" << src << " tar=" << tar << " len=" << len;	\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,	\
			  "Doing the memory transfer",	\
			  error, params.str());				\
    }									\
  }

#define CUDA_SET(tar, val, len, st)					\
  {									\
    cudaError_t error;							\
    std::stringstream params;						\
    if ((error = cudaMemsetAsync(tar, val, len, st)) != OK_CUDA){	\
      params << " tar=" << tar << " len=" << len;			\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,	\
			  "Doing the memory transfer",	\
			  error, params.str());				\
    }									\
  }


#define CUDA_FREE(ptr)							\
  {									\
    cudaError_t error;							\
    std::stringstream params;						\
    if ((error = cudaFree(ptr)) != OK_CUDA){				\
      params << "ptr=" << ptr;						\
      throw cudaException(__BASE_FILE__, __FUNCTION__, __LINE__,		\
			  "When the memory release", error, params.str()); \
    }									\
  }


#if TYPE_ID == FLOAT_ID
#define cusparseXcsrmv cusparseScsrmv
#define cublasXgemv cublasSgemv
#define cublasXot cublasSdot
#define cublasXaxpy cublasSaxpy
#define cublasXnorm cublasSnrm2

#elif TYPE_ID == DOUBLE_ID
#define cusparseXcsrmv cusparseDcsrmv
#define cublasXgemv cublasDgemv
#define cublasXot cublasDdot
#define cublasXaxpy cublasDaxpy
#define cublasXnorm cublasDnrm2

#else 
#error "The variable TYPE_ID is not defined"
#endif




#if DEBUG >= DEBUG_KERNEL_LVL
 #define KERNEL(func, dimG, dimB, stream, fmt, args...)		\
  _DEBUG(" Executing the kernel <%s> con:", #func);		\
  _DEBUG(" -Block mesh: [%d,%d,%d] ==> Total: %d",	\
	 dimG.x, dimG.y, dimG.z, dimG.x * dimG.y * dimG.z);	\
  _DEBUG(" -Threads mesh: [%d,%d,%d] ==> Total: %d",	\
	 dimB.x, dimB.y, dimB.z, dimB.x * dimB.y * dimB.z);	\
  _DEBUG(" -Param: "#fmt, args);				\
  func <<<dimG, dimB, 0, stream>>> (args);			\
  checkError(#func, __BASE_FILE__, __LINE__)
 #define KERNEL1(func, t1, dimG, dimB, stream, fmt, args...)	\
  _DEBUG("Executing the kernel <%s> con:", #func);		\
  _DEBUG(" -Block mesh: [%d,%d,%d] ==> Total: %d",	\
	 dimG.x, dimG.y, dimG.z, dimG.x * dimG.y * dimG.z);	\
  _DEBUG(" -Threads mesh: [%d,%d,%d] ==> Total: %d",	\
	 dimB.x, dimB.y, dimB.z, dimB.x * dimB.y * dimB.z);	\
  _DEBUG(" -Template: %s=%d", #t1, t1);				\
  _DEBUG(" -Param: "#fmt, args);				\
  func <t1> <<<dimG, dimB, 0, stream>>> (args);			\
  checkError(#func, __BASE_FILE__, __LINE__)
#else
 #define KERNEL(func, dimG, dimB, stream, fmt, args...)	\
  func <<<dimG, dimB, 0, stream>>> (args)

 #define KERNEL1(func, t1, dimG, dimB, stream, fmt, args...)	\
  func <t1> <<<dimG, dimB, 0, stream>>> (args)
#endif





#define TEXT_CUBLAS(x) (x==CUBLAS_STATUS_SUCCESS)? "OK":		\
  (x==CUBLAS_STATUS_NOT_INITIALIZED)? "NOT INITIALIZED":		\
  (x==CUBLAS_STATUS_ALLOC_FAILED)? "ALLOC FAILED":			\
  (x==CUBLAS_STATUS_INVALID_VALUE)? "INVALID VALUE":			\
  (x==CUBLAS_STATUS_ARCH_MISMATCH)? "ARCH MISMATCH":			\
  (x==CUBLAS_STATUS_MAPPING_ERROR)? "MAPPING FAILED":			\
  (x==CUBLAS_STATUS_EXECUTION_FAILED)? "EXECUTION FAILED":		\
  (x==CUBLAS_STATUS_INTERNAL_ERROR)? "INTERNAL ERROR":"OTHER"


#define TEXT_CUSPARSE(x) (x==CUSPARSE_STATUS_SUCCESS)? "OK":		\
  (x==CUSPARSE_STATUS_NOT_INITIALIZED)? "NOT INITIALIZED":		\
  (x==CUSPARSE_STATUS_ALLOC_FAILED)? "ALLOC FAILED":			\
  (x==CUSPARSE_STATUS_INVALID_VALUE)? "INVALID VALUE":			\
  (x==CUSPARSE_STATUS_ARCH_MISMATCH)? "ARCH MISMATCH":			\
  (x==CUSPARSE_STATUS_MAPPING_ERROR)? "MAPPING FAILED":			\
  (x==CUSPARSE_STATUS_EXECUTION_FAILED)? "EXECUTION FAILED":		\
  (x==CUSPARSE_STATUS_INTERNAL_ERROR)? "INTERNAL ERROR":		\
  (x==CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED)? "MATRIX TYPE NOT SUPPORTED":"OTHER"



#define TEXT_CUDA(x) (x==cudaSuccess)? "OK":				\
  (x==cudaErrorMissingConfiguration)?"cudaErrorMissingConfiguration" :	\
  (x==cudaErrorMemoryAllocation)?"cudaErrorMemoryAllocation" :		\
  (x==cudaErrorInitializationError)?"cudaErrorInitializationError" :	\
  (x==cudaErrorLaunchFailure)?"cudaErrorLaunchFailure" :		\
  (x==cudaErrorPriorLaunchFailure)?"cudaErrorPriorLaunchFailure" :	\
  (x==cudaErrorLaunchTimeout)?"cudaErrorLaunchTimeout" :		\
  (x==cudaErrorLaunchOutOfResources)?"cudaErrorLaunchOutOfResourc" :	\
  (x==cudaErrorInvalidDeviceFunction)?"cudaErrorInvalidDeviceFunction" : \
  (x==cudaErrorUnmapBufferObjectFailed)?"cudaErrorUnmapBufferObjectFailed" : \
  (x==cudaErrorInvalidDevice)?"cudaErrorInvalidDevice" :		\
  (x==cudaErrorInvalidValue)?"cudaErrorInvalidValue" :			\
  (x==cudaErrorInvalidPitchValue)?"cudaErrorInvalidPitchValue" :	\
  (x==cudaErrorInvalidSymbol)?"cudaErrorInvalidSymbol" :		\
  (x==cudaErrorMapBufferObjectFailed)?"cudaErrorMapBufferObjectFailed" : \
  (x==cudaErrorUnmapBufferObjectFailed)?"cudaErrorUnmapBufferObjectFailed" : \
  (x==cudaErrorInvalidHostPointer)?"cudaErrorInvalidHostPointer" :	\
  (x==cudaErrorInvalidDevicePointer)?"cudaErrorInvalidDevicePointer" :	\
  (x==cudaErrorInvalidTexture)?"cudaErrorInvalidTexture" :		\
  (x==cudaErrorInvalidTextureBinding)?"cudaErrorInvalidTextureBinding" : \
  (x==cudaErrorInvalidChannelDescriptor)?"cudaErrorInvalidChannelDescriptor" : \
  (x==cudaErrorInvalidMemcpyDirection)?"cudaErrorInvalidMemcpyDirection" : \
  (x==cudaErrorAddressOfConstant)?"cudaErrorAddressOfConstant" :	\
  (x==cudaErrorTextureFetchFailed)?"cudaErrorTextureFetchFailed" :	\
  (x==cudaErrorTextureNotBound)?"cudaErrorTextureNotBound" :		\
  (x==cudaErrorSynchronizationError)?"cudaErrorSynchronizationError" :	\
  (x==cudaErrorInvalidFilterSetting)?"cudaErrorInvalidFilterSetting" :	\
  (x==cudaErrorInvalidNormSetting)?"cudaErrorInvalidNormSetting" :	\
  (x==cudaErrorMixedDeviceExecution)?"cudaErrorMixedDeviceExecution" :	\
  (x==cudaErrorCudartUnloading)?"cudaErrorCudartUnloading" :		\
  (x==cudaErrorUnknown)?"cudaErrorUnknown" :				\
  (x==cudaErrorNotYetImplemented)?"cudaErrorNotYetImplemented" :	\
  (x==cudaErrorMemoryValueTooLarge)? "cudaErrorMemoryValueTooLarge":	\
  (x==cudaErrorInvalidResourceHandle)?"cudaErrorInvalidResourceHandle":	\
  (x==cudaErrorNotReady)?"cudaErrorNotReady":				\
  (x==cudaErrorInsufficientDriver)?"cudaErrorInsufficientDriver" :	\
  (x==cudaErrorSetOnActiveProcess)?"cudaErrorSetOnActiveProcess" :	\
  (x==cudaErrorInvalidSurface)?"cudaErrorInvalidSurface" :		\
  (x==cudaErrorNoDevice)?"cudaErrorNoDevice" :				\
  (x==cudaErrorECCUncorrectable)?"cudaErrorECCUncorrectable" :		\
  (x==cudaErrorSharedObjectSymbolNotFound)?"cudaErrorSharedObjectSymbolNotFound" : \
  (x==cudaErrorSharedObjectInitFailed)?"cudaErrorSharedObjectInitFailed" : \
  (x==cudaErrorUnsupportedLimit)?"cudaErrorUnsupportedLimit" :		\
  (x==cudaErrorDuplicateVariableName)?"cudaErrorDuplicateVariableName" : \
  (x==cudaErrorDuplicateTextureName)?"cudaErrorDuplicateTextureName" :	\
  (x==cudaErrorDuplicateSurfaceName)? "cudaErrorDuplicateSurfaceName":	\
  (x==cudaErrorDevicesUnavailable)?"cudaErrorDevicesUnavailable" :	\
  (x==cudaErrorInvalidKernelImage)?"cudaErrorInvalidKernelImage" :	\
  (x==cudaErrorNoKernelImageForDevice)?"cudaErrorNoKernelImageForDevice" : \
  (x==cudaErrorIncompatibleDriverContext)?"cudaErrorIncompatibleDriverContext" : \
  (x==cudaErrorPeerAccessAlreadyEnabled)?"cudaErrorPeerAccessAlreadyEnabled" : \
  (x==cudaErrorPeerAccessNotEnabled)?"cudaErrorPeerAccessNotEnabled" :	\
  (x==cudaErrorDeviceAlreadyInUse)?"cudaErrorDeviceAlreadyInUse" :	\
  (x==cudaErrorProfilerDisabled)?"cudaErrorProfilerDisabled" :		\
  (x==cudaErrorProfilerNotInitialized)?"cudaErrorProfilerNotInitialized" : \
  (x==cudaErrorProfilerAlreadyStarted)?"cudaErrorProfilerAlreadyStarted" : \
  (x==cudaErrorProfilerAlreadyStopped)?"cudaErrorProfilerAlreadyStopped" : \
  (x==cudaErrorStartupFailure)?"cudaErrorStartupFailure" :		\
  (x==cudaErrorApiFailureBase)?"cudaErrorApiFailureBase" : "other"



class cudaException : public std::exception {
 private:
  std::string msg;
  std::string type;
  std::string desc;
  std::string param;
  char * file;
  const char * func;
  int line;
  
 public:
  cudaException(char *file, const char *func, int line,
		const std::string& msg, cudaError_t error,
		const std::string& params);
  cudaException(char *file, const char *func, int line,
		const std::string& msg, cublasStatus_t error,
		const std::string& params);
  cudaException(char *file, const char *func, int line,
		const std::string& msg, cusparseStatus_t error,
		const std::string& params);
  virtual const char* what() const throw();
  void print() throw();
  std::string get_type() throw();
  std::string get_desc() throw();
  std::string get_param() throw();
  virtual ~cudaException() throw();
};


extern void checkError(const char *func, char *file, const int line);
#endif
