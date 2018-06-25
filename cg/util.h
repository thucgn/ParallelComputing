/*************************************************************************
	> File Name: util.h
	> Author: CGN
	> Mail: 
	> Created Time: 五  4/13 17:19:15 2018
 ************************************************************************/

#ifndef _UTIL_H
#define _UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>

#define LOG(format, ...) printf("[%s, %s:%d %s] " format "\n", \
                        __TIME__, __FILE__, __LINE__, __FUNCTION__, ##__VA_ARGS__)

#define LOGINFO(msg) printf("[%s, %s:%d %s] %s\n", \
                        __TIME__, __FILE__, __LINE__, __FUNCTION__, msg)

#define CHECK(flag) if(!(flag)) { \
                        printf("[%s, %s:%d %s] check failure!!!\n",\
                        __TIME__, __FILE__, __LINE__, __FUNCTION__);\
                        exit(0);}

#define MARK_TIME(t) gettimeofday(&t, NULL)
#define DIFF_TIME(tt, ts) (((tt).tv_sec-(ts).tv_sec) + ((tt).tv_usec*1e-6 - (ts).tv_usec*1e-6))

typedef struct timeval TIME_T;

#endif
