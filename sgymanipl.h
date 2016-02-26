/*
 * Frequency-domain Cadzow Fiter, Used to suppress random noise
 * of pre-stack seismic Data
 *
 * Copyright (C) 2013, Xi'an Jiaotong University Wenchao Chen
 * Create at: Fri 13 Sep 2013 12:10:34 PM CST 
 * Maintained by Wei Wang(email: weiwang.geo@gmail.com) 
 */

#ifndef __SEGYMANIPL_H__
#define __SEGYMANIPL_H__

#include "stdlib.h"
#include <stdio.h>
#include <string.h>
#include <math.h>


#define VOLUMEHEADLEN (3600)
#define TRACEHEADLEN (240)

struct _Segy_Volume_Header
{
	int Job;                                   //3204
	int Line;                                  //3208
	int Reel;                                  //3212

	short DataTracePerEnsemble;                //3214

	unsigned short AuxiliaryTracePerEnsemble;  //3216
	unsigned short dt;                         //3218
	unsigned short dtOrig;                     //3220
	unsigned short ns;                         //3222
	unsigned short nsOrig;                     //3224

	short DataSampleFormat;                    //3226
	short EnsembleFold;                        //3228
	short TraceSorting;                        //3230
	short VerticalSumCode;                     //3232

}; /* segy volume header */


struct _Segy_Trace_Header
{
	int TraceSequenceLine;             //4
	int TraceSequenceFile;             //8
	int Inline3D;                      //12
	int TraceNumber;                   //16
	int EnergySourcePoint;             //20
	int Crossline3D;                   //24
	int cdpTrace;                      //28

	int offset;                        //40

	int SourceX;                       //76
	int SourceY;                       //80
	int ReceiverX;                     //84
	int ReceiverY;                     //88

	short CoordinateUnits;             //90
	short WeatheringVelocity;          //92
	short SubWeatheringVelocity;       //94
	short SourceUpholeTime;            //96
	short GroupUpholeTime;             //98
	short SourceStaticCorrection;      //100
	short GroupStaticCorrection;       //102
	short TotalStaticApplied;          //104
	short LagTimeA;                    //106
	short LagTimeB;                    //108
	short DelayRecordingTime;          //110
	short MuteTimeStart;               //112
	short MuteTimeEND;                 //114

	unsigned short ns;                 //116
	unsigned short dt;                 //118

	int cdpX;                          //184
	int cdpY;                          //188
	int cdpElev;                       //192
	int FinalDatum;                    //196
	int ReplaceVel;                    //200
	int cdpStatics;                    //204
	int azimuth;                       //208
	int ovtNo;                         //212
	int nmoMute;                       //216

}; /* segy trace header */


int _get_volume_header(_Segy_Volume_Header &volume_header, unsigned char *volume_header_bytes);
int _put_volume_header(unsigned char *volume_header_bytes, _Segy_Volume_Header &volume_header);
int _get_trace_header(_Segy_Trace_Header &trace_header, unsigned char *trace_header_bytes);
int _put_trace_header(unsigned char *trace_header_bytes, _Segy_Trace_Header &trace_header);

int _get_bigendian_int(unsigned char *in_bytes, int *out_nums, int len);
int _get_bigendian_short(unsigned char *in_bytes, short *out_nums, int len);
int _get_bigendian_ushort(unsigned char *in_bytes, unsigned short *out_nums, int len);
int _put_bigendian_int(int *in_nums, unsigned char *out_bytes, int len);
int _put_bigendian_short(short *in_nums, unsigned char *out_bytes, int len);
int _put_bigendian_ushort(unsigned short *in_nums, unsigned char *out_bytes, int len);
int _IBM_BIG2IEEE_Float( unsigned char *in, float *out, int num);
int _IEEE_Float2IBM_BIG(float *input, unsigned char *output, int num);

#endif /* once __SGYMANIPL_H__ */