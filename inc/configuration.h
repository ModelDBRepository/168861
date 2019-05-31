//
//
// File author(s):  <Gabriel Rucabado and Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <time.h>
#include <string.h>
#include <stdio.h>
         
#define RANGE_ERROR 1E-14  
#define DOUBLE_ID 1
#define FLOAT_ID  2
 // Type of data [DOUBLE_ID || FLOAT_ID]    
#define TYPE_ID  1
#if (TYPE_ID == DOUBLE_ID)
#define TYPE double
#define TYPE_PRINT "%e"
#define strtoTYPE strtod
#elif (TYPE_ID == FLOAT_ID)
#define TYPE float
#define TYPE_PRINT "%e"
#define strtoTYPE strtod
#else
#error "Type of data does not exists"
#endif

#define SET_SEED srand48
#define GET_RAND drand48 
#define SIZE_NAME_LIB 1020
#define SCALE 0.001
#define RED "\033[22;31m"
#define BLACK "\033[22;30m"
#define GREEN "\033[22;32m"
#define BROWN "\033[22;33m"
#define BLUE "\033[22;34m"
#define MAGENTA "\033[22;35m"
#define CYAN "\033[22;36m"
#define GRAY "\033[22;37m"
#define YELLOW "\033[01;33m"
#define WHITE "\033[01;37m"
#define D_GRAY "\033[01;30m"
#define L_RED "\033[01;31m"
#define L_GREEN "\033[01;32m"
#define L_BLUE "\033[01;34m"
#define L_MAGENTA "\033[01;35m"
#define L_CYAN "\033[01;36m"
#define DEFAULT "\033[0m"

#define DEBUG_COLOR L_BLUE
#define INFO_COLOR GREEN
#define ERROR_COLOR RED


#define LOG(base, x) log(x)/log(base)
#define ABS(x) (((x)<0)?-(x):(x))
#define MIN(x, y) (((x)<(y))?(x):(y))


#define DEBUG_MESS_LVL 1

#if DEBUG == DEBUG_MESS_LVL
#define _DEBUG(fmt, args...)						\
  printf(DEBUG_COLOR"DEBUG: (%s) [%s -.- %s:%05d] ::: "fmt DEFAULT"\n",	\
	 Timestamp(),__BASE_FILE__, __FUNCTION__, __LINE__, ##args); fflush(stdout) 
#else
#define _DEBUG(d,f...)
#endif

#define _INFO_REP(fmt, args...)						\
  printf("\r"INFO_COLOR"INFO : (%s) ::: "fmt DEFAULT,			\
	 Timestamp(),##args); fflush(stdout) 
#define _INFO(fmt, args...)						\
  printf(INFO_COLOR"INFO : (%s) ::: "fmt DEFAULT"\n",			\
	 Timestamp(),##args); fflush(stdout) 
#define _ERROR(fmt, args...)						\
  fprintf(stderr,ERROR_COLOR"ERROR: (%s) [%s -.- %s:%05d] ::: "fmt DEFAULT"\n", \
	  Timestamp(),__BASE_FILE__, __FUNCTION__, __LINE__,##args); fflush(stdout)
  // For cuda
#define _ERROR2(file,func,line,fmt, args...)				\
  fprintf(stderr,ERROR_COLOR"ERROR: (%s) [%s -.- %s:%05d] ::: "fmt DEFAULT"\n",	\
	  Timestamp(), file, func, line,##args); fflush(stdout)
extern char *Timestamp();
#endif
