/*
 * Frequency-domain Cadzow Fiter, Used to suppress random noise
 * of pre-stack seismic Data
 *
 * Copyright (C) 2013, Xi'an Jiaotong University Wenchao Chen
 * Create at: Fri 13 Sep 2013 12:10:34 PM CST 
 * Maintained by Kai Yu (email: yukai.xjtu@gmail.com) 
 */

#ifndef __CADZOWFILTER_H__
#define __CADZOWFILTER_H__

#include "sgymanipl.h"
#include <stdio.h>
#include <memory.h>
#include <direct.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <complex>

//max pathname len
#define MAX_PATH 512
#define FFTMIN 512

using std::complex;
using std::size_t;

extern long _inlinePos, _xxlinePos;

//max pathname len
#define MAX_PATH 512
#define PI 3.1415926
#define EPS 1e-5

#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)


struct _dataPara
{
	bool is_big_endian;
	long inlineInt;
	long inlineEnd;
	long xxlineInt;
	long xxlineEnd;
	long tlen;
	double dt;
	long bkInsize;
	double maxAmplitude;
	double minAmplitude;
}; //data config




int _fio_stream(char *infilename, char *dipfilename, char *azimuthfilename, char *outfilename, const struct _dataPara &dpara);
void freqdomweights(char *input, char *output, const struct _dataPara &dpara);
void dirderiv2order(float *input, float *output, float *dip, float *azimuth, long *bkLen);
int _get_range(long *read_loc, long *write_loc, long max_len, long win_num, long win_len, long win_border);
size_t _get_fft_len(size_t len);
int _regular3d_format(char *ipathname, char *opathtrached, char *opathtracdat, _dataPara &dpara);
int _irregular3d_format(char *opathname, char *ipathtrached, char *ipathtracdat, const _dataPara &dpara);
void ArrayMaxMin(float arry[], int len, float *maxp, float *minp);

#endif /* once __CADZOWFILTER_H__ */

