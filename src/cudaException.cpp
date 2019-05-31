//
//
// File author(s):  <Gabriel Rucabado and Antonio Garcia Dopico>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file cudaExecption.cpp
  \brief Auxiliary function for the GPU implementation .
*/

#include <string.h>
#include "cudaException.h"
#include "configuration.h"
#include <exception>
#include <string.h>
#include <stdio.h>
#include <sstream>

extern void checkError(const char *func, char *file, const int line){
  cudaError_t error2;
  std::stringstream params;
    
  if ((error2 = cudaPeekAtLastError()) != OK_CUDA){
    params << "kernel=" << func << " row=" << file <<" line=" << line;
    throw cudaException(file, func, line,
			"Running kernel", error2, params.str());	
  }
}


const char* cudaException::what() const throw(){return NULL;}

cudaException::cudaException(char *file, const char *func, int line,
			     const std::string& msg, cublasStatus_t error,
			     const std::string& params){
  this->msg  	= msg;
  this->param	= params;
  this->type 	= TEXT_CUBLAS(error);
  this->desc 	= "Error cublas";
  this->file	= file;
  this->func	= func;
  this->line	= line;  
}
cudaException::cudaException(char *file, const char *func, int line,
			     const std::string& msg, cusparseStatus_t error,
			     const std::string& params){
  this->msg  	= msg;
  this->param	= params;
  this->type 	= TEXT_CUSPARSE(error);
  this->desc 	= "Error cusparse";
  this->file	= file;
  this->func	= func;
  this->line	= line;
}
cudaException::cudaException(char *file, const char *func, int line,
			     const std::string& msg, cudaError_t error,
			     const std::string& params) {
  this->msg  	= msg;
  this->param	= params;
  this->type 	= TEXT_CUDA(error);
  this->desc 	= cudaGetErrorString(error);
  this->file	= file;
  this->func	= func;
  this->line	= line;
 }

cudaException::~cudaException() throw(){}
void cudaException::print() throw(){
  _ERROR2(file, func, line, "Captured exception:");
  _ERROR2(file, func, line, 
	 " -Type: %s"
	 " -Desc: %s"
	 " -Param: %s",
	 type.c_str(), desc.c_str(), param.c_str());
}
std::string cudaException::get_type() throw(){return type;}
std::string cudaException::get_desc() throw(){return desc;}
std::string cudaException::get_param() throw(){return param;}
