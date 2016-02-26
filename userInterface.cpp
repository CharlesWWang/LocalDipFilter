/****
* Copyright (c) 2015, University of Texas at Austin
* All rights reserved.
* 
* Version: 1.0
* Authors: Wei Wang, Yuhong Guo
* Date   : 12/28/2015
****/

#include "localdipfilter.h"
#include "sgymanipl.h"

long _inlinePos, _xxlinePos;

int main(int argc, char *argv[])
{
	// 定义输入输出文件指针
    FILE *fpinf;
	// 定义输入输出文件名
	char infilename[MAX_PATH];
	char outfilename[MAX_PATH];
	char indfilename[MAX_PATH];
	char inafilename[MAX_PATH];
	char inffilename[MAX_PATH];
	errno_t ferr, err;
	int temp;


	//define filter & data parameters
	struct _dataPara dpara;

	fprintf(stdout, "***********LOCAL DIP FILTER USING PRECALCULATED DIP&AZIMUTH!***********\n");
	fprintf(stdout, "The program is designed to calculate apparent dips for 3D seismic data.\n");
	fprintf(stdout, "The input data should be regular, all traces required with same length.\n");
	fprintf(stdout, "File 'dipcalc3dinformation.txt' need be created in curr project folder.\n\n");
	
	//get information file pointer
	getcwd(inffilename, MAX_PATH);
	strcat(inffilename, "\\localdipfilterinformation.txt");
	fprintf(stdout, "Information file name is:\n %s.\n", inffilename);
	ferr = fopen_s(&fpinf, inffilename, "r");	
    if (ferr){
		fprintf(stderr, "Information file read error: %d! Check or Create if not exist!\n", ferr);
        goto prret;
    }

	//read input file name including path
	char tempstr[MAX_PATH];
	fgets(tempstr, MAX_PATH ,fpinf);
	if (!feof(fpinf)){
		fgets(infilename, MAX_PATH, fpinf);
		infilename[strlen(infilename) - 1] = '\0';
		fprintf(stdout, "Input data    file name is:\n %s\n", infilename);
	}

	fgets(tempstr, MAX_PATH ,fpinf);
	if (!feof(fpinf)){
		fgets(indfilename, MAX_PATH, fpinf);
		indfilename[strlen(indfilename) - 1] = '\0';
		fprintf(stdout, "Input dip     file name is:\n %s\n", indfilename);
	}

	fgets(tempstr, MAX_PATH ,fpinf);
	if (!feof(fpinf)){
		fgets(inafilename, MAX_PATH, fpinf);
		inafilename[strlen(inafilename) - 1] = '\0';
		fprintf(stdout, "Input azimuth file name is:\n %s\n", inafilename);
	}

	fgets(tempstr, MAX_PATH ,fpinf);
	if (!feof(fpinf)){
		fgets(outfilename, MAX_PATH, fpinf);
		outfilename[strlen(outfilename) - 1] = '\0';
		fprintf(stdout, "Output data   file name is:\n %s\n", outfilename);
	}
	
	//read file format
	if (!feof(fpinf)){
		fgets(tempstr, MAX_PATH, fpinf);	
	}
	if(!feof(fpinf)){
		fscanf(fpinf, "%s\n", tempstr);
		temp = atof(tempstr);
		switch(temp) 
		{
			case 0:
				dpara.is_big_endian = false;
				break;
			case 1:
				dpara.is_big_endian = true;
				break;
			default:
				dpara.is_big_endian = true;
		}
		fprintf(stderr, "The datfile format is: %ld\n", temp);
    }
		
	//read inline range
	if (!feof(fpinf)){
		fgets(tempstr, MAX_PATH, fpinf);	
	}
	if(!feof(fpinf)){
		fscanf(fpinf, "%s\n", tempstr);
		dpara.inlineInt = atof(tempstr);
		fprintf(stdout, "The initial inline is: %ld\n", dpara.inlineInt);
		fscanf(fpinf, "%s\n", tempstr);
		dpara.inlineEnd = atof(tempstr);
		fprintf(stdout, "The last inline is   : %ld\n", dpara.inlineEnd);	
	}

	//read xxline range
	if (!feof(fpinf)){
		fgets(tempstr, MAX_PATH, fpinf);	
	}
	if(!feof(fpinf)){
		fscanf(fpinf, "%s\n", tempstr);
		dpara.xxlineInt = atof(tempstr);
		fprintf(stdout, "The initial xxline is: %ld\n", dpara.xxlineInt);
		fscanf(fpinf, "%s\n", tempstr);
		dpara.xxlineEnd = atof(tempstr);
		fprintf(stdout, "The last xxline is   : %ld\n", dpara.xxlineEnd);	
	}

	//read time length & sampling rate
	if (!feof(fpinf)){
		fgets(tempstr, MAX_PATH, fpinf);	
	}
	if(!feof(fpinf)){
		fscanf(fpinf, "%s\n", tempstr);
		dpara.tlen = atof(tempstr);
		fprintf(stdout, "The time length is   : %ld\n", dpara.tlen);
		fscanf(fpinf, "%s\n", tempstr);
		dpara.dt = atof(tempstr);
		fprintf(stdout, "The sampling rate is : %lf s\n", dpara.dt);	
	}

	//read inline position in trace header
	if (!feof(fpinf)){
		fgets(tempstr, MAX_PATH, fpinf);	
	}
	if(!feof(fpinf)){
		fscanf(fpinf, "%s\n", tempstr);
		_inlinePos = atof(tempstr);
		fprintf(stdout, "The inline is restored in trace header: %ld\n", _inlinePos);
	}

	//read xxline position in trace header
	if (!feof(fpinf)){
		fgets(tempstr, MAX_PATH, fpinf);	
	}
	if(!feof(fpinf)){
		fscanf(fpinf, "%s\n", tempstr);
		_xxlinePos = atof(tempstr);
		fprintf(stdout, "The xxline is restored in trace header: %ld\n", _xxlinePos);
	}

	//read dip calculation parameters subset inline numbers
	if (!feof(fpinf)){
		fgets(tempstr, MAX_PATH, fpinf);	
	}
	if(!feof(fpinf)){
		fscanf(fpinf, "%s\n", tempstr);
		dpara.bkInsize = atof(tempstr);
	}
	long maxInlineNum = dpara.inlineEnd - dpara.inlineInt  + 1;
	if (dpara.bkInsize>maxInlineNum)
		dpara.bkInsize = maxInlineNum;
	else if (dpara.bkInsize<=0){  
		fprintf(stderr, "Block size inline dim should be positive bigger than filter size!\n");
		goto prret;
	}
	fprintf(stdout, "Subblocks for Inline subset process is : %ld\n", dpara.bkInsize);


	//define interprocess temporary data
	char intracedat[MAX_PATH];
	char intracehed[MAX_PATH];
	char outracedat[MAX_PATH];
	char indiphed[MAX_PATH];
	char indip[MAX_PATH];
	char inazimuth[MAX_PATH];
	getcwd(intracedat, MAX_PATH);
	getcwd(intracehed, MAX_PATH);
	getcwd(outracedat, MAX_PATH);
	getcwd(indiphed, MAX_PATH);
	getcwd(indip, MAX_PATH);
	getcwd(inazimuth, MAX_PATH);
	strcat(intracedat, "\\LocalDipFiltOrig.dat");
	strcat(intracehed, "\\LocalDipFiltHead.dat");
	strcat(outracedat, "\\LocalDipFiltSmth.dat");
	strcat(indiphed, "\\LocalDipFiltDIPHead.dat");
	strcat(indip, "\\LocalDipFiltDIPP.dat");
	strcat(inazimuth, "\\LocalDipFiltAZIM.dat");

	//transform file format
	fprintf(stdout, "\nLoading original SEGY data ...\n");
	ferr = _regular3d_format(infilename, intracehed, intracedat, dpara);
    if (ferr){
		fprintf(stderr, "Original 3d volume regular format error: %d! Check if exist!\n", ferr);
        goto prret;
    }

	//transform file format
	fprintf(stdout, "\nLoading dip SEGY data ...\n");
	ferr = _regular3d_format(indfilename, indiphed, indip, dpara);
    if (ferr){
		fprintf(stderr, "Dip 3d volume regular format error: %d! Check if exist!\n", ferr);
        goto prret;
    }

	//transform file format
	fprintf(stdout, "\nLoading azimuth SEGY data ...\n");
	ferr = _regular3d_format(inafilename, indiphed, inazimuth, dpara);
    if (ferr){
		fprintf(stderr, "Azimuth 3d volume regular format error: %d! Check if exist!\n", ferr);
        goto prret;
    }

	//原始数据体分块的块数
	fprintf(stdout, "\n******LOCAL DIP FILTERING START******\n");

	//开始分块处理原始数据
	//_fio_stream(intracedat, indip, inazimuth, outracedat, dpara);

	
	//transform file format(careful!! Output now stored in INTRACEDAT)
	ferr = _irregular3d_format(outfilename, outracedat, intracehed, dpara);
	if (ferr){
		fprintf(stderr, "Filtering Segy file restoration error: %d! Check if it exist!\n", ferr);
        goto prret;
    }
/*
	err = remove(intracehed);
	if(err){
		fprintf(stderr, "Temporary File (Original Trace Head) Deleting Error!\n");
		goto prret;
	}

	err = remove(intracedat);
	if(err){
		fprintf(stderr, "Temporary File (Original Trace Data) Deleting Error!\n");
		goto prret;
	}

	err = remove(indip);
	if(err){
		fprintf(stderr, "Temporary File (Dip Data) Deleting Error!\n");
		goto prret;
	}

	err = remove(inazimuth);
	if(err){
		fprintf(stderr, "Temporary File (Azimuth Data) Deleting Error!\n");
		goto prret;
	}

	err = remove(outracedat);
	if(err){
		fprintf(stderr, "Temporary File (Filtered Data) Deleting Error!\n");
		goto prret;
	}
	*/

	// 完成原始数据的滤波处理
prret:
	fprintf(stdout, "\n****************************************************\n");
	fprintf(stdout, "******DIPS CALCULATION USING WVDF METHOD FINISH!******\n");
	fprintf(stdout, "******************************************************\n");
    system("pause");
}

