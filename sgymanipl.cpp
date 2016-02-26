/*
 * This file forms the main part of segy file manipulations
 *
 * Copyright (C) 2015, Xi'an Jiaotong University
 * Maintained by Wei Wang (email: weiwang.geo@gmail.com) 
 */

#include "sgymanipl.h"
#include "localdipfilter.h"

/*
_get_volume_header: read segy volume header information into the continuous bytes string of 
                    the format of unsigned char in big endiant  
*/
int _get_volume_header(_Segy_Volume_Header &volume_header, unsigned char *volume_header_bytes)
{
	int chunkInt[3];
	short chunkShort[4];
	unsigned short chunkUShort[5];

	_get_bigendian_int(&volume_header_bytes[3200], chunkInt, 3);          //3212
	volume_header.Job = chunkInt[0];
	volume_header.Line = chunkInt[1];
	volume_header.Reel = chunkInt[2];

	_get_bigendian_short(&volume_header_bytes[3212], chunkShort, 1);      //3214
	volume_header.DataTracePerEnsemble = chunkShort[0];

	_get_bigendian_ushort(&volume_header_bytes[3214], chunkUShort, 5);    //3224
	volume_header.AuxiliaryTracePerEnsemble = chunkUShort[0];
	volume_header.dt = chunkUShort[1];
	volume_header.dtOrig = chunkUShort[2];
	volume_header.ns = chunkUShort[3];
	volume_header.nsOrig = chunkUShort[4];

	_get_bigendian_short(&volume_header_bytes[3224], chunkShort, 4);      //3232
	volume_header.DataSampleFormat = chunkShort[0];
	volume_header.EnsembleFold = chunkShort[1];
	volume_header.TraceSorting = chunkShort[2];
	volume_header.VerticalSumCode = chunkShort[3];

	return 0;
}

/*
_put_volume_header: write segy volume header information into the continuous bytes string of 
                   the format of unsigned char in big endiant  
*/
int _put_volume_header(unsigned char *volume_header_bytes, _Segy_Volume_Header &volume_header)
{
	int chunkInt[3];
	short chunkShort[4];
	unsigned short chunkUShort[5];

	chunkInt[0] = volume_header.Job;
	chunkInt[1] = volume_header.Line;
	chunkInt[2] = volume_header.Reel;
	_put_bigendian_int(chunkInt, &volume_header_bytes[3200], 3);

	chunkShort[0] = volume_header.DataTracePerEnsemble;
	_put_bigendian_short(chunkShort, &volume_header_bytes[3212], 1);

	chunkUShort[0] = volume_header.AuxiliaryTracePerEnsemble;
	chunkUShort[1] = volume_header.dt;
	chunkUShort[2] = volume_header.dtOrig;
	chunkUShort[3] = volume_header.ns;
	chunkUShort[4] = volume_header.nsOrig;
	_put_bigendian_ushort(chunkUShort, &volume_header_bytes[3214], 5);

	chunkShort[0] = volume_header.DataSampleFormat;
	chunkShort[1] = volume_header.EnsembleFold;
	chunkShort[2] = volume_header.TraceSorting;
	chunkShort[3] = volume_header.VerticalSumCode;
	_put_bigendian_short(&volume_header.DataSampleFormat, &volume_header_bytes[3224], 4);

	return 0;
}

/*
_get_trace_header: read nessary segy trace header information from the continuous bytes string  
                   of the format of unsigned char in big endiant  
*/
int _get_trace_header(_Segy_Trace_Header &trace_header, unsigned char *trace_header_bytes)
{
	int chunkInt[9];
	short chunkShort[13];
	unsigned short chunkUShort[2];

	_get_bigendian_int(trace_header_bytes, chunkInt, 7);                 //28
	trace_header.TraceSequenceLine = chunkInt[0];
	trace_header.TraceSequenceFile = chunkInt[1];
	trace_header.Inline3D = chunkInt[2];
	trace_header.TraceNumber = chunkInt[3];
	trace_header.EnergySourcePoint = chunkInt[4];
	trace_header.Crossline3D = chunkInt[5];
	trace_header.cdpTrace = chunkInt[6];                                     

	_get_bigendian_int(&trace_header_bytes[36], chunkInt, 1);            //40
	trace_header.offset = chunkInt[0];

	_get_bigendian_int(&trace_header_bytes[72], chunkInt, 4);
	trace_header.SourceX = chunkInt[0];
	trace_header.SourceY = chunkInt[1];
	trace_header.ReceiverX = chunkInt[2];
	trace_header.ReceiverY = chunkInt[3];                                //88

	_get_bigendian_short(&trace_header_bytes[88], chunkShort, 13);       //114
	trace_header.CoordinateUnits = chunkShort[0];
	trace_header.WeatheringVelocity = chunkShort[1];
	trace_header.SubWeatheringVelocity = chunkShort[2];
	trace_header.SourceUpholeTime = chunkShort[3];
	trace_header.GroupUpholeTime = chunkShort[4];
	trace_header.SourceStaticCorrection = chunkShort[5];
	trace_header.GroupStaticCorrection = chunkShort[6];
	trace_header.TotalStaticApplied = chunkShort[7];
	trace_header.LagTimeA = chunkShort[8];
	trace_header.LagTimeB = chunkShort[9];
	trace_header.DelayRecordingTime = chunkShort[10];
	trace_header.MuteTimeStart = chunkShort[11];
	trace_header.MuteTimeEND = chunkShort[12];                            

	_get_bigendian_ushort(&trace_header_bytes[114], chunkUShort, 2);
	trace_header.ns = chunkUShort[0];
	trace_header.dt = chunkUShort[1];                                    //118

	_get_bigendian_int(&trace_header_bytes[180], chunkInt, 9);           //216
	trace_header.cdpX = chunkInt[0];
	trace_header.cdpY = chunkInt[1];
	trace_header.cdpElev = chunkInt[2];
	trace_header.FinalDatum = chunkInt[3];
	trace_header.ReplaceVel = chunkInt[4];
	trace_header.cdpStatics = chunkInt[5];
	trace_header.azimuth = chunkInt[6];
	trace_header.ovtNo = chunkInt[7];
	trace_header.nmoMute = chunkInt[8];  

	//for special trace header puting
	_get_bigendian_int(&trace_header_bytes[_inlinePos-1], chunkInt, 1);
	trace_header.Inline3D = chunkInt[0];
	_get_bigendian_int(&trace_header_bytes[_xxlinePos-1], chunkInt, 1);
	trace_header.Crossline3D = chunkInt[0];

	return 0;
}

/*
_put_trace_header: write segy trace header information into the continuous bytes string of 
                   the format of unsigned char in big endiant  
*/
int _put_trace_header(unsigned char *trace_header_bytes, _Segy_Trace_Header &trace_header)
{
	int chunkInt[9];
	short chunkShort[13];
	unsigned short chunkUShort[2];
	
	chunkInt[0] = trace_header.TraceSequenceLine;
	chunkInt[1] = trace_header.TraceSequenceFile;
	chunkInt[2] = trace_header.Inline3D;
	chunkInt[3] = trace_header.TraceNumber;
	chunkInt[4] = trace_header.EnergySourcePoint;
	chunkInt[5] = trace_header.Crossline3D;
	chunkInt[6] = trace_header.cdpTrace;                                     
	_put_bigendian_int(chunkInt, trace_header_bytes, 7);               //28

	chunkInt[0] = trace_header.offset;
	_put_bigendian_int(chunkInt, &trace_header_bytes[36], 1);          //40

	chunkInt[0] = trace_header.SourceX;
	chunkInt[1] = trace_header.SourceY;
	chunkInt[2] = trace_header.ReceiverX;
	chunkInt[3] = trace_header.ReceiverY;                                 
	_put_bigendian_int(chunkInt, &trace_header_bytes[72], 4);          //88

	chunkShort[0] = trace_header.CoordinateUnits;
	chunkShort[1] = trace_header.WeatheringVelocity;
	chunkShort[2] = trace_header.SubWeatheringVelocity;
	chunkShort[3] = trace_header.SourceUpholeTime;
	chunkShort[4] = trace_header.GroupUpholeTime;
	chunkShort[5] = trace_header.SourceStaticCorrection;
	chunkShort[6] = trace_header.GroupStaticCorrection;
	chunkShort[7] = trace_header.TotalStaticApplied;
	chunkShort[8] = trace_header.LagTimeA;
	chunkShort[9] = trace_header.LagTimeB;
	chunkShort[10] = trace_header.DelayRecordingTime;
	chunkShort[11] = trace_header.MuteTimeStart;
	chunkShort[12] = trace_header.MuteTimeEND;                            
	_put_bigendian_short(chunkShort, &trace_header_bytes[88], 13);     //114

	chunkUShort[0] = trace_header.ns;
	chunkUShort[1] = trace_header.dt;                                     
	_put_bigendian_ushort(chunkUShort, &trace_header_bytes[114], 2);   //118

	chunkInt[0] = trace_header.cdpX;
	chunkInt[1] = trace_header.cdpY;
	chunkInt[2] = trace_header.cdpElev;
	chunkInt[3] = trace_header.FinalDatum;
	chunkInt[4] = trace_header.ReplaceVel;
	chunkInt[5] = trace_header.cdpStatics;
	chunkInt[6] = trace_header.azimuth;
	chunkInt[7] = trace_header.ovtNo;
	chunkInt[8] = trace_header.nmoMute;                                   
	_put_bigendian_int(chunkInt, &trace_header_bytes[180], 9);         //216


	//for special trace header puting
	chunkInt[0] = trace_header.Inline3D;
	_put_bigendian_int(chunkInt, &trace_header_bytes[_inlinePos-1], 1);
	chunkInt[0] = trace_header.Crossline3D;
	_put_bigendian_int(chunkInt, &trace_header_bytes[_xxlinePos-1], 1);
	
	return 0;
}

/*
_get_bigendian_int: transform INT numbers of big endiant in the continuous bytes string of 
                    unsigned char format into little endiant  
*/
int _get_bigendian_int(unsigned char *in_bytes, int *out_nums, int len)
{
	int i; 
	double sign, number;
	unsigned char byte1, byte2, byte3, byte4;

	for (i = 0; i < len; ++i){
		byte1 = *in_bytes++;
		byte2 = *in_bytes++;
		byte3 = *in_bytes++;
		byte4 = *in_bytes++;

		sign = (double)(byte1 >> 7); // gain sign from first bit
		number = (double)(byte1&0x7f);
		number = (((number * 256.0) + byte2) * 256.0 + byte3) * 256.0 + byte4;

		out_nums[i] = (int)(1 - 2 * sign) * number; 	
	}

	return 0;
}

/*
_get_bigendian_short: transform SHORT numbers of big endiant in the continuous bytes string of 
                      unsigned char format into little endiant  
*/
int _get_bigendian_short(unsigned char *in_bytes, short *out_nums, int len)
{
	int i;
	double sign, number;
	unsigned char byte1, byte2;

	for (i = 0; i < len; ++i){
		byte1 = *in_bytes++;
		byte2 = *in_bytes++;

		sign = (double)(byte1 >> 7);
		number = (double)(byte1&0x7f);
		number = (number * 256.0) + byte2;
		out_nums[i] = (short)(1 - 2 * sign) * number;	
	}

	return 0;
}

/*
_get_bigendian_ushort: transform UNSIGNED SHORT numbers of big endiant in the continuous bytes 
                       string of unsigned char format into little endiant  
*/
int _get_bigendian_ushort(unsigned char *in_bytes, unsigned short *out_nums, int len)
{
	int i;
	unsigned char byte1, byte2;

	for (i = 0; i < len; ++i){
		byte1 = *in_bytes++;
		byte2 = *in_bytes++;

		out_nums[i] = (unsigned short)((byte1 * 256.0) + byte2);	
	}

	return 0;
}

/*
_put_bigendian_int: transform INT numbers of little endiant into the continuous bytes 
                    string of unsigned char format  
*/
int _put_bigendian_int(int *in_nums, unsigned char *out_bytes, int len)
{
	int i;
	int number, sign;

	for (i = 0; i < len; ++i){
		number = in_nums[i];
		sign = ( number < 0 ? 1 : 0);
		number = number * (1 - 2 * sign); // abs(input)
		if (number > 0){ //nonzeros
			out_bytes[4 * i + 3] = (unsigned char)(number % 256);
			number = number / 256;
			out_bytes[4 * i + 2] = (unsigned char)(number % 256);
			number = number / 256;
			out_bytes[4 * i + 1] = (unsigned char)(number % 256);
			out_bytes[4 * i] = (unsigned char)(number/256 + sign * 128);
		}
		else
		{
			out_bytes[4 * i] = (unsigned char)0;
			out_bytes[4 * i + 1] = (unsigned char)0;		
			out_bytes[4 * i + 2] = (unsigned char)0;		
			out_bytes[4 * i + 3] = (unsigned char)0;
			continue;
		}
	}

	return 0;
}

/*
_put_bigendian_short: transform SHORT numbers of little endiant into the continuous bytes 
                      string of unsigned char format  
*/
int _put_bigendian_short(short *in_nums, unsigned char *out_bytes, int len)
{
	int i;
	short number, sign;

	for (i = 0; i < len; ++i){
		number = in_nums[i];
		sign = ( number < 0 ? 1 : 0);
		number = number * (1 - 2 * sign); // abs(input)
		if (number > 0){ //nonzeros
			out_bytes[2 * i + 1] = (unsigned char)(number % 256);
			number = number / 256;
			out_bytes[2 * i] = (unsigned char)(number + sign * 128);
		}
		else
		{
			out_bytes[2 * i] = (unsigned char)0;
			out_bytes[2 * i + 1] = (unsigned char)0;		
			continue;
		}
	}

	return 0;
}

/*
_put_bigendian_ushort: transform UNSIGNED SHORT numbers of little endiant into the continuous bytes 
                       string of unsigned char format  
*/
int _put_bigendian_ushort(unsigned short *in_nums, unsigned char *out_bytes, int len)
{
	int i;
	short number;

	for (i = 0; i < len; ++i){
		number = in_nums[i];
		if (number > 0){ //nonzeros
			out_bytes[2 * i + 1] = (unsigned char)(number % 256);
			out_bytes[2 * i] = (unsigned char)(number / 256);
		}
		else
		{
			out_bytes[2 * i] = (unsigned char)0;
			out_bytes[2 * i + 1] = (unsigned char)0;		
			continue;
		}
	}

	return 0;
}

/*
_IBM_BIG2IEEE_Float: transform signle-precision floating-point of big endiant in forms of the continuous bytes 
                       string of unsigned char into little endiant IEEE format  
*/
int _IBM_BIG2IEEE_Float( unsigned char *in, float *out, int num)
{
    int m;
    
    unsigned char byte1,byte2,byte3,byte4;
    
    double sign, expo, frac;
    double denominator = (double)pow(2.0,24);
    for( m = 0; m<num; m++ )
    {
        byte1 = *in++;
        byte2 = *in++; 
        byte3 = *in++;    
        byte4 = *in++;
        
        sign = (double)( byte1 >>7) ; // gain sign from first bit 
        expo = (double)( (byte1&0x7f) ) - 64;// gain exponent from first byte,
        frac = ( ( byte2*256 + byte3 )*256 + byte4 )/denominator; // gain mantissa from last 3 bytes 
        
        out[m] = (float)( 1-2*sign)*( pow(16 ,expo) ) * frac;
     } 
   
    return 0;
         
}

/*
_IEEE_Float2IBM_BIG: transform signle-precision floating-point of little endiant in forms of the continuous bytes 
                       string of unsigned char into big endiant IBM format  
*/
int _IEEE_Float2IBM_BIG(float *input, unsigned char *output, int num)
{
    int m;
	long sign;   //sign 
	long exp;    //exponent
	long fmant;  //fraction
    float datain, tmp ; // attention : cannot use   long input1;
	double fm;
     
	for(m=0;m<num;m++) 
	{ 
		datain = input[m];
		sign = ( datain<0?1:0 ) ;
		datain = datain*(1-2*sign);// abs(input)
		exp = 0;
		tmp = datain;

		if (datain>0)    // nonzeros
		{
			if(datain>=1.0)
			{
				exp++;
				while( tmp/16.f >= 1.0)
				{
					exp++;
					tmp= tmp/16.f;
				}
			}  
			else
			{
				while ( tmp*16.f < 1.0 )
				{
					exp--;
					tmp=tmp*16.f;
				}
				//exp++;// attention : ibm fmant   :     0.mant   not 1.mant !
			}

		}
		else
		{
			exp = 0;
			fmant = (long)0;
			output[4*m] = (unsigned char)exp;
			output[4*m+3] = (unsigned char)0;		
			output[4*m+2] = (unsigned char)0;		
			output[4*m+1] = (unsigned char)0;
			continue;
		}
    
		exp = ( exp + 64 );                 //exponent part 
		fm = datain * pow(16.0,-(exp-64));  //fraction of floating-point
		fmant=(long)( fm * pow(2.0,24) ) ;  //尾数无符号整型
		//dataout = (sign<<31)|(exp<<24)|fmant ;  
		output[4*m] = (unsigned char)exp + sign*128;
		output[4*m+3] = (unsigned char)fmant%256;
		fmant = ( fmant - output[4*m+3] )/256;
		output[4*m+2] = (unsigned char)fmant%256;
		fmant = ( fmant - output[4*m+2] )/256;
		output[4*m+1] = (unsigned char)fmant;
	}	
	return 0;
}


