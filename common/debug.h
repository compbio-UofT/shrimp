/*
 * debug.h
 *
 *  Created on: 2009-11-06
 *      Author: dlister
 */

#ifndef DEBUG_H_
#define DEBUG_H_
//uncomment to debug
#define DEBUG(x,...) fprintf(stderr,"DEBUG: " x "\n",##__VA_ARGS__)
#ifndef DEBUG
#define DEBUG(X,...)
#endif

#endif /* DEBUG_H_ */
