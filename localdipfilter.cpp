/*
 * Frequency-domain Cadzow Fiter, Used to suppress random noise
 * of pre-stack seismic Data
 *
 * Copyright (C) 2013, Xi'an Jiaotong University Wenchao Chen
 * Create at: Fri 13 Sep 2013 12:10:34 PM CST 
 * Maintained by Kai Yu (email: yukai.xjtu@gmail.com) 
 */

#include "localdipfilter.h"
#pragma comment(lib, "libfftw3-3.lib")


//delete the fftw plan
#define del_plan(x) do { if (x) { fftw_destroy_plan(x); (x) = NULL; } } while(0)


//**************************************************************
//_FIO_STREAM, constructed to prepare one proper size dataset for
//             filtering, store processed data to disk file
//**************************************************************
int _fio_stream(char *infilename, char *dipfilename, char *azimuthfilename, char *outfilename, const struct _dataPara &dpara)
{
	long double_byte = sizeof(double);
	long float_byte = sizeof(float);
	long char_byte = sizeof(char);
	errno_t err = 0;
	time_t start_time, finish_time;

	//determine filtering direction and virtual data size
	long maxInlineNum, maxXxlineNum;
	maxInlineNum = dpara.inlineEnd - dpara.inlineInt + 1; 
	maxXxlineNum = dpara.xxlineEnd - dpara.xxlineInt + 1; 

	//prepare file names including path
	char preDiff[MAX_PATH];
	getcwd(preDiff, MAX_PATH);
	strcat(preDiff, "\\temporary_pre_diff.dat");
    //perform frequency domain coordinates weighting
	freqdomweights(infilename, preDiff, dpara);
	//////////////////////////////////////////////////////////////////////////

	long lenY, lenX, lenZ, bkLenY, bkLenX, bkLenZ;
	long gborder, overlap, ny;
	long bkLen[3];
	lenY = maxInlineNum;
	lenX = maxXxlineNum;
	lenZ = dpara.tlen;
	bkLenY = dpara.bkInsize;
	bkLenX = lenX;
	bkLenZ = lenZ;
	bkLen[0] = bkLenY;
	bkLen[1] = bkLenX;
	bkLen[2] = bkLenZ;

	gborder = 1;
	overlap = 2 * gborder;
	//calculate no of block division
	long bkNumY = ceil((float)(maxInlineNum-overlap)/(bkLenY-overlap));
	fprintf(stdout, "Block division number along Inline is: %ld\n", bkNumY);

	//计算每个分块沿inline方向读取和存储位置(free() followed)
	long *readInlineLoc = new long[2*bkNumY+1];
	long *writInlineLoc = new long[2*bkNumY+1];
	_get_range(readInlineLoc, writInlineLoc, maxInlineNum, bkNumY, bkLenY, gborder);

	float *orig_data = NULL;
	float *diff_data = NULL;
	float *dip = NULL;
	float *azimuth = NULL;
	orig_data = (float*) malloc(bkLenY * bkLenX * bkLenZ * float_byte);
	diff_data = (float*) malloc(bkLenY * bkLenX * bkLenZ * float_byte);
	dip = (float*) malloc(bkLenY * bkLenX * bkLenZ * float_byte);
	azimuth = (float*) malloc(bkLenY * bkLenX * bkLenZ * float_byte);
	FILE *ifp, *ofp, *ifp_dip, *ifp_azimuth;
    if (fopen_s(&ofp, outfilename, "wb")){
		fprintf(stderr, "Fatal error! Output filtering file open error!\n");
        goto fret;
    }
    if (fopen_s(&ifp, preDiff, "rb")){
		fprintf(stderr, "Fatal error! Output pre diff file open error!\n");
        goto fret;
    }
    if (fopen_s(&ifp_dip, dipfilename, "rb")){
		fprintf(stderr, "Fatal error! Output filtering file open error!\n");
        goto fret;
    }
    if (fopen_s(&ifp_azimuth, azimuthfilename, "rb")){
		fprintf(stderr, "Fatal error! Output filtering file open error!\n");
        goto fret;
    }

	long xx, yy, zz;
	__int64 lineOffset, dataOffset;
	for (ny=0; ny<bkNumY; ++ny){
		//将文件指针移动到当前分块数据的起始侧线
		lineOffset = readInlineLoc[ny]*lenX*lenZ*float_byte;
		if(_fseeki64(ifp, lineOffset, SEEK_SET)){
			perror("Fatal Error: fseek failed!");
			break;
		}
		if(_fseeki64(ifp_dip, lineOffset, SEEK_SET)){
			perror("Fatal Error: fseek failed!");
			break;
		}
		if(_fseeki64(ifp_azimuth, lineOffset, SEEK_SET)){
			perror("Fatal Error: fseek failed!");
			break;
		}
		//初始化数据指针
		start_time = time(NULL);
		dataOffset = 0;
		for (yy=readInlineLoc[ny]; yy<=readInlineLoc[ny+bkNumY]; yy++)
			for (xx=0; xx<lenX; xx++){
				fread((orig_data+dataOffset),float_byte,lenZ,ifp);
				fread((dip+dataOffset),float_byte,lenZ,ifp_dip);
				fread((azimuth+dataOffset),float_byte,lenZ,ifp_azimuth);
				dataOffset = dataOffset + lenZ;
			}

		dirderiv2order(orig_data, diff_data, dip, azimuth, bkLen);

		//将文件指针移动到当前分块数据起始侧线
		lineOffset = writInlineLoc[ny]*lenX*lenZ*float_byte;
		if(_fseeki64(ofp, lineOffset, SEEK_SET)){
			perror("Fatal Error: fseek failed!");
			break;
		}
		for (yy=writInlineLoc[ny]; yy<=writInlineLoc[ny+bkNumY]; yy++){
			dataOffset = (yy-readInlineLoc[ny])*lenX*lenZ;
			for (xx=0; xx<lenX; xx++){
				fwrite(diff_data+dataOffset,float_byte,lenZ,ofp);
				dataOffset += lenZ;
			}
		}
		finish_time = time(NULL);
		fprintf(stdout, "%ld / %ld FINISHED USING TIME %lf s!\n", 
				(ny+1), bkNumY, difftime(finish_time, start_time));
	}

fret:
	free(orig_data);
	free(diff_data);
	free(dip);
	free(azimuth);

	fclose(ifp);
	fclose(ifp_dip);
	fclose(ifp_azimuth);
	fclose(ofp);

	/*
	err = remove(preDiff);
	if(err){
		fprintf(stderr, "Temporary File (Frequency Coordinates Weighting Data) Deleting Error!\n");
		goto fret;
	}
	*/

    return 0;
}


/*******************************************************************
freqdomweights: Frequency Domain Coordinates Weighting
	input - input data
	output - output coordinates weighted time domain data
	dpara - struct containing information about data
********************************************************************/
void freqdomweights(char *input, char *output, const struct _dataPara &dpara)
{
	long double_byte = sizeof(double);
	long float_byte = sizeof(float);
	long char_byte = sizeof(char);
	errno_t err = 0;
	
	//timer
	time_t start_time, finish_time;
	time_t wstart_time, wfinish_time;

	wstart_time = time(NULL);

	//define fftw plans
	fftw_plan forward_plan_1d = NULL;
	fftw_plan reverse_plan_1d = NULL;
	fftw_plan forward_plan_2d = NULL;
	fftw_plan reverse_plan_2d = NULL;

	//define 2d & 1d fft2 stacks   
	complex<double> *fftw_freq = NULL;
	complex<double> *fftw_time = NULL;
	complex<double> *fftw_time_1d = NULL;
	complex<double> *fftw_freq_1d = NULL;

	//define data stacks
	float *trace_buff = NULL;
	double *temp_fft_real = NULL;
	double *temp_fft_imag = NULL;

	/////////////////////////////////////////////////////////////
	//determine filtering direction and virtual data size
	long maxInlineNum, maxXxlineNum;
	maxInlineNum = dpara.inlineEnd - dpara.inlineInt + 1; 
	maxXxlineNum = dpara.xxlineEnd - dpara.xxlineInt + 1; 
	//get the closest power of 2
	long fft_lenz, fft_lenx, fft_leny, fft_len_max, half_fft_lenz, half_fft_len_max;
	fft_lenz = _get_fft_len(dpara.tlen);
	half_fft_lenz = fft_lenz / 2 + 1;
	fft_lenx = _get_fft_len(maxXxlineNum);
	fft_leny = _get_fft_len(maxInlineNum);
	fft_len_max = MAX(fft_lenx, fft_leny);
	half_fft_len_max = fft_lenz / 2 + 1;

	//assign stack space
	fftw_time = (complex<double>*) fftw_malloc(sizeof(complex<double>) * fft_lenz * fft_len_max);
	fftw_freq = (complex<double>*) fftw_malloc(sizeof(complex<double>) * fft_lenz * fft_len_max);
	fftw_time_1d = (complex<double>*) fftw_malloc(sizeof(complex<double>) * fft_len_max);
	fftw_freq_1d = (complex<double>*) fftw_malloc(sizeof(complex<double>) * fft_len_max);
	trace_buff = (float*) malloc(float_byte * dpara.tlen);//float
	temp_fft_real = (double*) malloc(double_byte * fft_lenz);
	temp_fft_imag = (double*) malloc(double_byte * fft_lenz);

	//intialize fftw plan
    forward_plan_1d = fftw_plan_dft_1d(fft_len_max, reinterpret_cast<fftw_complex*>(fftw_time_1d),
            reinterpret_cast<fftw_complex*>(fftw_freq_1d), FFTW_FORWARD, FFTW_MEASURE);
	reverse_plan_1d = fftw_plan_dft_1d(fft_len_max, reinterpret_cast<fftw_complex*>(fftw_freq_1d),
            reinterpret_cast<fftw_complex*>(fftw_time_1d), FFTW_BACKWARD, FFTW_MEASURE);

	forward_plan_2d = fftw_plan_dft_2d(fft_lenz, fft_len_max, reinterpret_cast<fftw_complex*>(fftw_time),
		    reinterpret_cast<fftw_complex*>(fftw_freq), FFTW_FORWARD, FFTW_MEASURE);
	reverse_plan_2d = fftw_plan_dft_2d(fft_lenz, fft_len_max, reinterpret_cast<fftw_complex*>(fftw_freq), 
		    reinterpret_cast<fftw_complex*>(fftw_time), FFTW_BACKWARD, FFTW_MEASURE);

	//prepare file names including path
	char ifftData[MAX_PATH];
	char fftReal[MAX_PATH];
	char fftImag[MAX_PATH];
	getcwd(ifftData, MAX_PATH);
	getcwd(fftReal, MAX_PATH);
	getcwd(fftImag, MAX_PATH);
	strcat(ifftData, "\\temporary_ifft.dat");
	strcat(fftReal, "\\temporary_fft_real.dat");
	strcat(fftImag, "\\temporary_fft_imag.dat");
    //open file stream to start processing
	FILE *ifp, *ofp, *ofp_real, *ofp_imag;
    if (fopen_s(&ifp, input, "rb")){
		fprintf(stderr, "Fatal error! Input original file open error!\n");
        goto wret;
    }
	//store intermediary data pre-2nd order difference with float format
    if (fopen_s(&ofp, output, "wb+")){
		fprintf(stderr, "Fatal error! Output pre diff file open error!\n");
        goto wret;
    }
	//store intermediary complex data of fft, each part is double format
    if (fopen_s(&ofp_real, fftReal, "wb+")){
		fprintf(stderr, "Fatal error! Output real part file open error!\n");
        goto wret;
    }	
    if (fopen_s(&ofp_imag, fftImag, "wb+")){
		fprintf(stderr, "Fatal error! Output imag part file open error!\n");
        goto wret;
    }

	start_time = time(NULL);
	long xx, yy, zz;
	memset(fftw_time, 0, sizeof(complex<double>) * fft_lenz * fft_len_max);
	//2d fft applied along xxline&time direction
	for (yy=0; yy<maxInlineNum; ++yy){
		for (xx=0; xx<maxXxlineNum; ++xx)
			if (fread(trace_buff, float_byte, dpara.tlen, ifp)==dpara.tlen){
				for (zz=0; zz<dpara.tlen; ++zz)
					fftw_time[xx*fft_lenz+zz] = complex<double>(trace_buff[zz], 0.0);
			}else{
				fprintf(stderr, "Fatal error! Input original file read error!\n");
			}
		
		//2d fft
		fftw_execute_dft(forward_plan_2d, reinterpret_cast<fftw_complex*>(fftw_time), reinterpret_cast<fftw_complex*>(fftw_freq));

		//temporarily write results to disk
		//size - maxInlineNum * fft_len_max * fft_lenz
		for (xx=0; xx<fft_len_max; ++xx){
			for (zz=0; zz<fft_lenz; ++zz){
				temp_fft_real[zz] = fftw_freq[xx*fft_lenz+zz].real();
				temp_fft_imag[zz] = fftw_freq[xx*fft_lenz+zz].imag();
			}
			fwrite(temp_fft_real, double_byte, fft_lenz, ofp_real);
			fwrite(temp_fft_imag, double_byte, fft_lenz, ofp_imag);
		}
	}
	finish_time = time(NULL);
	fprintf(stdout, "2d-FFTs implemented along Xxline are finished for all Inlines!\n");
	fprintf(stdout, "THE TOTALLY COST TIME IS %lf s!\n", difftime(finish_time, start_time));

	//1d fft + 2d ifft applied along inline direction
	start_time = time(NULL);
	long fft_trace_len_bytes = fft_lenz * double_byte;
	__int64 rtrace_beg_offset, rtrace_end_offset;
	double norm;
	memset(fftw_time_1d, 0, sizeof(complex<double>) * fft_len_max);
	if(_fseeki64(ofp_real, (__int64)0, SEEK_SET))
		fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
	if(_fseeki64(ofp_imag, (__int64)0, SEEK_SET))
		fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
	for (yy=0; yy<fft_len_max; ++yy){
		if(_fseeki64(ofp_real, (__int64)0, SEEK_SET))
			fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
		if(_fseeki64(ofp_imag, (__int64)0, SEEK_SET))
			fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
		rtrace_beg_offset = yy * fft_trace_len_bytes;
		rtrace_end_offset = (fft_len_max-1-yy) * fft_trace_len_bytes;
		for (xx=0; xx<maxInlineNum; ++xx){	
			if(_fseeki64(ofp_real, rtrace_beg_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
			if(_fseeki64(ofp_imag, rtrace_beg_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
			if(fread(temp_fft_real, double_byte, fft_lenz, ofp_real) != fft_lenz)
				fprintf(stderr, "Fatal error! ofp_real reading failed!\n");
			if(fread(temp_fft_imag, double_byte, fft_lenz, ofp_imag) != fft_lenz)
				fprintf(stderr, "Fatal error! ofp_imag reading failed!\n");
			for (zz=0; zz<fft_lenz; ++zz)
				fftw_time[xx*fft_lenz+zz] = complex<double>(temp_fft_real[zz], temp_fft_imag[zz]);	
			if(_fseeki64(ofp_real, rtrace_end_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
			if(_fseeki64(ofp_imag, rtrace_end_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
		}

		for (xx=0; xx<fft_lenz; ++xx){
			for (zz=0; zz<maxInlineNum; ++zz)
				fftw_time_1d[zz] = fftw_time[zz*fft_lenz+xx];
			//1d fft
			fftw_execute_dft(forward_plan_1d, reinterpret_cast<fftw_complex*>(fftw_time_1d), reinterpret_cast<fftw_complex*>(fftw_freq_1d));
			for (zz=0; zz<fft_len_max; ++zz)
				fftw_freq[zz*fft_lenz+xx] = fftw_freq_1d[zz];
		}

		//**************************************************************************
		//index norm weighting
		norm = (2*PI*yy/fft_len_max)*(2*PI*yy/fft_len_max);
		for (xx=1; xx<fft_len_max; ++xx){
			norm += (2*PI*xx/fft_len_max)*(2*PI*xx/fft_len_max);
			for (zz=1; zz<fft_lenz; ++zz){					
				norm += (2*PI*zz/fft_lenz)*(2*PI*zz/fft_lenz) + EPS;
				fftw_freq[xx*fft_len_max+zz] = fftw_freq[xx*fft_len_max+zz]/norm;
			}
		}
		//**************************************************************************

		//2d ifft
		fftw_execute_dft(reverse_plan_2d, reinterpret_cast<fftw_complex*>(fftw_freq), reinterpret_cast<fftw_complex*>(fftw_time));
		if(_fseeki64(ofp_real, (__int64)0, SEEK_SET))
			fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
		if(_fseeki64(ofp_imag, (__int64)0, SEEK_SET))
			fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
		rtrace_beg_offset = yy * fft_trace_len_bytes;
		rtrace_end_offset = (fft_len_max-1-yy) * fft_trace_len_bytes;
		for (xx=0; xx<maxInlineNum; ++xx){	
			if(_fseeki64(ofp_real, rtrace_beg_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
			if(_fseeki64(ofp_imag, rtrace_beg_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
			for (zz=0; zz<fft_lenz; ++zz){
				temp_fft_real[zz] = fftw_time[xx*fft_lenz+zz].real();
				temp_fft_imag[zz] = fftw_time[xx*fft_lenz+zz].imag();
			}
			fwrite(temp_fft_real, double_byte, fft_lenz, ofp_real);
			fwrite(temp_fft_imag, double_byte, fft_lenz, ofp_imag);
			if(_fseeki64(ofp_real, rtrace_end_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
			if(_fseeki64(ofp_imag, rtrace_end_offset, SEEK_CUR))
				fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
		}

		//fprintf(stdout, "2d-IFFTs implemented along Inline %ld are finished!\n", yy);
	}
	finish_time = time(NULL);
	fprintf(stdout, "1d-FFTs implemented along Inline are finished for all Xxlines!\n");
	fprintf(stdout, "2d-IFFTs implemented along Inline are finished for all Xxlines!\n");
	fprintf(stdout, "THE TOTALLY COST TIME IS %lf s!\n", difftime(finish_time, start_time));

	//1d ifft applied along xxline direction
	if(_fseeki64(ofp_real, (__int64)0, SEEK_SET))
		fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_real!\n");
	if(_fseeki64(ofp_imag, (__int64)0, SEEK_SET))
		fprintf(stderr, "Fatal error! _fseeki64 failed for ofp_imag!\n");
	start_time = time(NULL);
	for (yy=0; yy<maxInlineNum; ++yy){
		for (xx=0; xx<fft_len_max; ++xx){
			if (fread(temp_fft_real, double_byte, fft_lenz, ofp_real)!=fft_lenz)
				fprintf(stderr, "Fatal error! ofp_real reading failed!\n");
			if (fread(temp_fft_imag, double_byte, fft_lenz, ofp_imag)!=fft_lenz)
				fprintf(stderr, "Fatal error! ofp_imag reading failed!\n");
			for (zz=0; zz<fft_lenz; ++zz)
				fftw_time[xx*fft_lenz+zz] = complex<double>(temp_fft_real[zz], temp_fft_imag[zz]);
		}

		for (xx=0; xx<fft_lenz; ++xx){
			for (zz=0; zz<fft_len_max; ++zz)
				fftw_freq_1d[zz] = fftw_time[zz*fft_lenz+xx];
			//1d ifft
			fftw_execute_dft(reverse_plan_1d, reinterpret_cast<fftw_complex*>(fftw_freq_1d), reinterpret_cast<fftw_complex*>(fftw_time_1d));
			for (zz=0; zz<fft_len_max; ++zz)
				fftw_time[zz*fft_lenz+xx] = fftw_time_1d[zz];		
		}

		for (xx=0; xx<maxXxlineNum; ++xx){
			for (zz=0; zz<dpara.tlen; ++zz)
				trace_buff[zz] = fftw_time[xx*fft_lenz+zz].real()/fft_lenz/fft_len_max/fft_len_max;
			fwrite(trace_buff, float_byte, dpara.tlen, ofp);
		}
	}
	finish_time = time(NULL);
	fprintf(stdout, "1d-iFFTs implemented along Xxline are finished for all Inlines!\n");
	fprintf(stdout, "THE TOTALLY COST TIME IS %lf s!\n", difftime(finish_time, start_time));
	fprintf(stdout, "************3D FFTs and 3D IFFTs are all finished!*************\n");
	wfinish_time = time(NULL);
	fprintf(stdout, "THE TOTALLY COST TIME IS %lf s!\n", difftime(wfinish_time, wstart_time));

	fclose(ifp);
	fclose(ofp);
	fclose(ofp_real);
	fclose(ofp_imag);

	err = remove(fftReal);
	if(err){
		fprintf(stderr, "Temporary File (FFT Real Part) Deleting Error!\n");
		goto wret;
	}

	err = remove(fftImag);
	if(err){
		fprintf(stderr, "Temporary File (FFT Imag Part) Deleting Error!\n");
		goto wret;
	}
	
wret:
	del_plan(forward_plan_1d);
    del_plan(reverse_plan_1d);
	del_plan(forward_plan_2d);
    del_plan(reverse_plan_2d);

	fftw_free(fftw_time);
	fftw_free(fftw_freq);
	fftw_free(fftw_time_1d);
	fftw_free(fftw_freq_1d);

	free(trace_buff);
	free(temp_fft_real);
	free(temp_fft_imag);

	fclose(ifp);
	fclose(ofp_real);
	fclose(ofp_imag);
	fclose(ofp);
}


/*******************************************************************
isograd: Isotropic Gradient Operator
	input - input data
	dt, dx, dy - output gradient components
	bkLen - size of data (0-inline, 1-xxline, 2-time)
********************************************************************/
void dirderiv2order(float *input, float *output, float *dip, float *azimuth, long *bkLen)
{
	//get&put data pointer
	float *pdip, *pazimuth, *pin, *pout; 
	float nt, nx, ny, ddip, dazimuth;
	//border effect of iso-gradient
	long gborder = 1;
	long Nt, Nx, Ny, Nxt, N;
	Ny = bkLen[0]; 
	Nx = bkLen[1];
	Nt = bkLen[2];
	Nxt = Nx * Nt;
	N = Nxt * Ny;
	
	long yy, xx, tt;
	//initializing output
	pout = output;
	pin = input;
	for (yy=0; yy<Ny; ++yy)
		for (xx=0; xx<Nx; ++xx)
			for (tt=0; tt<Nt; ++tt){
				if (yy<gborder || xx<gborder || tt<gborder || yy>Ny-gborder-1 || xx>Nx-gborder-1 || tt>Nt-gborder-1)
					*pout = *pin;
				pout++;
				pin ++;
			}

	//move input pointer to skip borders
	pin = input + gborder + gborder * Nt + gborder * Nxt;
	pout = output + gborder + gborder * Nt + gborder * Nxt;
	pdip = dip + gborder + gborder * Nt + gborder * Nxt;
	pazimuth = azimuth + gborder + gborder * Nt + gborder * Nxt;
	for (yy = 0; yy < Ny-2*gborder; ++yy){
		for(xx = 0; xx <Nx-2*gborder; ++xx){
			for (tt = 0; tt <Nt-2*gborder; ++tt){
				ddip = (*pdip) * PI / 180;
				dazimuth = (*pazimuth) * PI / 180;
				nt = sin(ddip);
				nx = cos(ddip) * cos(dazimuth);
				ny = cos(ddip) * sin(dazimuth);
				//input[tt][xx][yy] = *(pin + tt + xx * Nt + yy * Nxt)
				*pout = (*pin) * (-2*nx*nx - 2*ny*ny - 2*nt*nt); //[yy][xx][tt]
				*pout += *(pin + Nt) * (nx * nx); //[yy][xx + 1][tt]
				*pout += *(pin - Nt) * (nx * nx); //[yy][xx - 1][tt]
				*pout += *(pin + Nxt) * (ny * ny); //[yy + 1][xx][tt]
				*pout += *(pin - Nxt) * (ny * ny); //[yy - 1][xx][tt]
				*pout += *(pin + 1) * (nt * nt); //[yy][xx][tt + 1]
				*pout += *(pin - 1) * (nt * nt); //[yy][xx][tt - 1]
				*pout += *(pin + Nt + Nxt) * (nx * ny / 2); //[yy + 1][xx + 1][tt]
				*pout += *(pin - Nt - Nxt) * (nx * ny / 2); //[yy - 1][xx - 1][tt]
				*pout += *(pin - Nt + Nxt) * (-nx * ny / 2); //[yy + 1][xx - 1][tt]
				*pout += *(pin + Nt - Nxt) * (-nx * ny / 2); //[yy - 1][xx + 1][tt]
				*pout += *(pin + Nt + 1) * (nx * nt / 2); //[yy][xx + 1][tt + 1]
				*pout += *(pin - Nt - 1) * (nx * nt / 2); //[yy][xx - 1][tt - 1]
				*pout += *(pin - Nt + 1) * (-nx * nt / 2); //[yy][xx - 1][tt + 1]
				*pout += *(pin + Nt - 1) * (-nx * nt / 2); //[yy][xx + 1][tt - 1]
				*pout += *(pin + Nxt + 1) * (ny * nt / 2); //[yy + 1][xx][tt + 1]
				*pout += *(pin - Nxt - 1) * (ny * nt / 2); //[yy - 1][xx][tt - 1]
				*pout += *(pin - Nxt + 1) * (-ny * nt / 2); //[yy - 1][xx][tt + 1]
				*pout += *(pin + Nxt - 1) * (-ny * nt / 2); //[yy + 1][xx][tt - 1]

				*pout *= (-1.0);

				//move to next analysis point
				pin ++;
				pout ++;
				pdip ++;
				pazimuth ++;				
			}
			//move to next xxline
			pin += 2 * gborder;
			pout += 2 * gborder;
			pdip += 2 * gborder;
			pazimuth += 2 * gborder;
		}
		//move to next inline
		pin  += 2 * gborder * Nt;
		pout += 2 * gborder * Nt;
		pdip += 2 * gborder * Nt;
		pazimuth += 2 * gborder * Nt;
	}
}


/****************************************************************
_get_range: calculate input & output positions when filtering with border effect
	read_loc  2*lx, for reading positions,
	               [0][...] store front loc, [1][...] store rear loc
	write_loc 2*lx, for writing positions,
	               [0][...] store front loc, [1][...] store rear loc
	max_len, the maximum size of data along some dimension
	win_num, the number of total window division
	win_len, the input window size of a filter
	win_border, the number of influenced points due to border effect
*****************************************************************/
int _get_range(long *read_loc, long *write_loc, long max_len, long win_num, long win_len, long win_border)
{
	long temp, win_exp_len, t;
	//the size to export in each window
	win_exp_len = win_len - 2 * win_border;
	temp = win_len;


	if (win_num == 1){	
		read_loc[0] = 0;
		read_loc[0 + win_num] = max_len - 1;
		write_loc[0] = 0;
		write_loc[0 + win_num] = max_len - 1;
		return 0;
	}
	else if (win_num > 1){
		for (t=0; t<win_num; ++t){//read&write the first subset
			if (t == 0){
				read_loc[t] = 0;
				read_loc[t+win_num] = win_len - 1;

				write_loc[t] = 0;
				write_loc[t+win_num] = win_len - 1 - win_border;			
			}
			else if (t == win_num - 1){//read&write the last subset
				read_loc[t] = max_len - win_len;
				read_loc[t+win_num] = max_len - 1;

				write_loc[t] = max_len - win_len + win_border;
				write_loc[t+win_num] = max_len - 1;
			}
			else{//read&write middle subsets
				temp = (win_len - 2 * win_border) * t;
				read_loc[t] = temp;
				read_loc[t+win_num] = temp + win_len - 1;
	
				write_loc[t] = temp + win_border;
				write_loc[t+win_num] = temp + win_len - 1 - win_border;	
			}
		}
		return 0;	
	}
	else{
		return -1;
	}
}


/****************************************************************
 * ensure the fft length is more than 1024
 * and must be power of 2
 ****************************************************************/
size_t _get_fft_len(size_t len)
{
    size_t i;

    if (len < FFTMIN)
        return FFTMIN;

    for (i = 1; i < len; i *= 2)
        ; 

    return i;
}

/***********************************************************************************
_REGULAR3D_FORMAT separate segy data into trace data and segy header, and trace data
are converted from IBM format to IEEE format, segy volume_header & trace_header are stored 
in one file sequentially, besides non-full inlines are added zero traces to get regular
3d volume.
************************************************************************************/
int _regular3d_format(char *ipathname, char *opathtrached, char *opathtracdat, _dataPara &dpara)
{
	long FloatBytes = sizeof(float);
	long CharBytes = sizeof(char);
	int err = 0;

    //open input segy file
	FILE *ifp;
	errno_t ferr = fopen_s(&ifp, ipathname, "rb");	
    if (ferr){
		fprintf(stderr, "Fatal error! Input file (_regular3d_format) read error: %d\n", ferr);
		err = -1;
        goto cfret;
    }
	//open output data file
	FILE *ofp_d;
	ferr = fopen_s(&ofp_d, opathtracdat, "wb");	
    if (ferr){
		fprintf(stderr, "Fatal error! Data file (_regular3d_format) output error: %d\n", ferr);
		fclose(ifp);
		err = -1;
        goto cfret;
    }
	//open output header file
	FILE *ofp_h;
	ferr = fopen_s(&ofp_h, opathtrached, "wb");
    if (ferr){
		fprintf(stderr, "Fatal error! Header file (_regular3d_format) output error: %d\n", ferr);
		fclose(ifp);
		fclose(ofp_d);
		err = -1;
        goto cfret;
    }

	//transfer segy volume header
	unsigned char VolumeHeadChars[VOLUMEHEADLEN];
	fread(VolumeHeadChars, CharBytes, VOLUMEHEADLEN, ifp);
	fwrite(VolumeHeadChars, CharBytes, VOLUMEHEADLEN, ofp_h);
	//transfer segy traces
	unsigned char TraceHeadChars[TRACEHEADLEN];
	unsigned char NullTraceHeadChars[TRACEHEADLEN];
	float *TraceBuff = NULL;	
	unsigned char *TraceBuffChars = NULL;
	long TraceLengthChars = FloatBytes * dpara.tlen;
	TraceBuff = (float*)malloc(TraceLengthChars);
	TraceBuffChars = (unsigned char*)malloc(TraceLengthChars);

	struct _Segy_Trace_Header TraceHeaderSt;
	struct _Segy_Trace_Header NullTraceHeaderSt;
	long WritInline, WritXxline, ReadInline, ReadXxline;
	dpara.maxAmplitude = 0;
	dpara.minAmplitude = 0;
	float maxTemp, minTemp;
	long readTag = 1;
	for (WritInline=dpara.inlineInt; WritInline<=dpara.inlineEnd; ++WritInline){
		for (WritXxline=dpara.xxlineInt; WritXxline<=dpara.xxlineEnd; ++WritXxline){
			if (readTag){//try to read one trace
				err = fread(TraceHeadChars, CharBytes, TRACEHEADLEN, ifp);
				if (err == TRACEHEADLEN){
					_get_trace_header(TraceHeaderSt, TraceHeadChars);						
					ReadInline = TraceHeaderSt.Inline3D;
					ReadXxline = TraceHeaderSt.Crossline3D;
				}else{
					fprintf(stderr, "Read Segy File Trace Header (Inline=%ld, Xxline=%ld) Fatal Error!\n", WritInline, WritXxline);
				}
				err = fread(TraceBuffChars, CharBytes, TraceLengthChars, ifp);
				if (err == TraceLengthChars){
					_IBM_BIG2IEEE_Float(TraceBuffChars, TraceBuff, dpara.tlen);
					ArrayMaxMin(TraceBuff, dpara.tlen, &maxTemp, &minTemp);
					dpara.maxAmplitude = dpara.maxAmplitude > maxTemp ? dpara.maxAmplitude : maxTemp;
					dpara.minAmplitude = dpara.minAmplitude < minTemp ? dpara.minAmplitude : minTemp;
				}else{
					fprintf(stderr, "Read Segy File Trace Data (Inline=%ld, Xxline=%ld) Fatal Error!\n", WritInline, WritXxline);
				}
			}
			if (WritInline==ReadInline && WritXxline==ReadXxline){
				fwrite(TraceHeadChars, CharBytes, TRACEHEADLEN, ofp_h);
				fwrite(TraceBuff, FloatBytes, dpara.tlen, ofp_d);
				readTag = 1;			
			}else{
				NullTraceHeaderSt = TraceHeaderSt;
				NullTraceHeaderSt.Inline3D = WritInline;
				NullTraceHeaderSt.Crossline3D = WritXxline;
				_put_trace_header(NullTraceHeadChars, NullTraceHeaderSt);
				fwrite(NullTraceHeadChars, CharBytes, TRACEHEADLEN, ofp_h);
				fwrite(TraceBuff, FloatBytes, dpara.tlen, ofp_d);
				fprintf(stdout, "Null Trace (Inline=%ld, Xxline=%ld) Is Filled by Trace (Inline=%ld, Xxline=%ld)!\n", WritInline, WritXxline, ReadInline, ReadXxline);
				readTag = 0;			
			}		
		}//xxline	
	}//inline

	err = 0; //case 0: all rightly

cfret:
	if (TraceBuffChars){
		free(TraceBuffChars);
	}
	if (TraceBuff){
		free(TraceBuff);
	}
	if (ifp){
		fclose(ifp);
	}
	if (ofp_h){
		fclose(ofp_h);
	}
	if (ofp_d){
		fclose(ofp_d);
	}
    return err;
}



/***********************************************************************************
_IRREGULAR3D_FORMAT compose trace data and segy header to obtain SUN format segy data,
trace data in IEEE format need to be converted into IBM format, segy volume_header &
trace_header are stored in one file.
************************************************************************************/
int _irregular3d_format(char *opathname, char *ipathtrached, char *ipathtracdat, const _dataPara &dpara)
{
	float *TraceBuff = NULL;	
	unsigned char *TraceBuffChars = NULL;

	long FloatBytes = sizeof(float);
	long CharBytes = sizeof(char);
	int err = 0;

    //open input segy file
	FILE *ofp;
	errno_t ferr = fopen_s(&ofp, opathname, "wb");	
    if (ferr){
		fprintf(stderr, "Fatal error! Iutput file read error: %d\n", ferr);
		err = -1;
        goto dfret;
    }
	//open input data file
	FILE *ifp_d;
	ferr = fopen_s(&ifp_d, ipathtracdat, "rb");	
    if (ferr){
		fprintf(stderr, "Fatal error! Output data file read error: %d\n", ferr);
		fclose(ofp);
		err = -1;
        goto dfret;
    }
	//open input header file
	FILE *ifp_h;
	ferr = fopen_s(&ifp_h, ipathtrached, "rb");
    if (ferr){
		fprintf(stderr, "Fatal error! Iutput file read error: %d\n", ferr);
		fclose(ofp);
		fclose(ifp_d);
		err = -1;
        goto dfret;
    }

	//transfer segy volume header
	unsigned char VolumeHeadChars[VOLUMEHEADLEN];
	fread(VolumeHeadChars, CharBytes, VOLUMEHEADLEN, ifp_h);
	fwrite(VolumeHeadChars, CharBytes, VOLUMEHEADLEN, ofp);
	//transfer segy traces
	unsigned char TraceHeadChars[TRACEHEADLEN];

	long TraceLengthChars = FloatBytes * dpara.tlen;
	TraceBuff = (float*)malloc(TraceLengthChars);
	TraceBuffChars = (unsigned char*)malloc(TraceLengthChars);
	struct _Segy_Trace_Header TraceHeaderSt;
	long WritInline, WritXxline, ReadInline, ReadXxline;
	for (WritInline=dpara.inlineInt; WritInline<=dpara.inlineEnd; ++WritInline){
		for (WritXxline=dpara.xxlineInt; WritXxline<=dpara.xxlineEnd; ++WritXxline){
			err = fread(TraceHeadChars, CharBytes, TRACEHEADLEN, ifp_h);
			if (err == TRACEHEADLEN){
				_get_trace_header(TraceHeaderSt, TraceHeadChars);						
				ReadInline = TraceHeaderSt.Inline3D;
				ReadXxline = TraceHeaderSt.Crossline3D;
			}else{
				fprintf(stderr, "Read Segy Trace Header File (Inline=%ld, Xxline=%ld) Fatal Error!\n", WritInline, WritXxline);
			}//read trace header of each trace
			err = fread(TraceBuff, FloatBytes, dpara.tlen, ifp_d);
			if (err == dpara.tlen){
				_IEEE_Float2IBM_BIG(TraceBuff, TraceBuffChars, dpara.tlen);
			}else{
				fprintf(stderr, "Read Segy Trace Data File (Inline=%ld, Xxline=%ld) Fatal Error!\n", WritInline, WritXxline);
			}//read trace data of each trace

			if (WritInline==ReadInline && WritXxline==ReadXxline){
				fwrite(TraceHeadChars, CharBytes, TRACEHEADLEN, ofp);
				fwrite(TraceBuffChars, CharBytes, TraceLengthChars, ofp);
			}else{
				err = 2;
				fprintf(stderr, "Segy Trace Header File Isn't Composed Properly!\n");
				fprintf(stdout,"(Inline=%ld, Xxline=%ld) Is Read from Segy Trace Header File!\n", ReadInline, ReadXxline);
				fprintf(stdout,"(Inline=%ld, Xxline=%ld) Should Be Writen!\n", WritInline, WritXxline);
				goto dfret;	
			}		
		}//xxline	
	}//inline

	err = 0; //case 0: all rightly

dfret:
	if (TraceBuffChars){
		free(TraceBuffChars);
	}
	if (TraceBuff){
		free(TraceBuff);
	}
	if (ofp){
		fclose(ofp);
	}
	if (ifp_h){
		fclose(ifp_h);
	}
	if (ifp_d){
		fclose(ifp_d);
	}
    return err;
}


/****************************************************************
ArrayMaxMin: Search for the minima and maxima of a float arry
*****************************************************************/
void ArrayMaxMin(float arry[], int len, float *maxp, float *minp)
{
	int i;
	float temp;
	*maxp = arry[0];
	*minp = *maxp;
	for(i=1; i<=len-1; i++)
	{
		temp = arry[i];
		if(temp >= *maxp)
			*maxp = temp;
		if(temp < *minp)
			*minp = temp;
	}
}