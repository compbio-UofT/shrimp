/*
 * debug.h
 *
 *  Created on: 2009-11-06
 *      Author: dlister
 */

#ifndef DEBUG_H_
#define DEBUG_H_

//#define DEBUGGING
//#define DEBUG_MAIN
//#define DEBUG_SEEDS

#ifdef DEBUGGING
#define DEBUG(x,...) fprintf(stderr,"DEBUG: " x "\n",##__VA_ARGS__)
#else
#define DEBUG(X,...)
#endif

#endif /* DEBUG_H_ */
