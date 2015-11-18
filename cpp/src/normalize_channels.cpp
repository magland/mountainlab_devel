#include "normalize_channels.h"
#include <stdio.h>
#include "mdaio.h"
#include <math.h>
#include <QTime>

void normalize_channels_1(MDAIO_HEADER &H,double *sum,double *sumsqr,int NN,FILE *input_file);
void normalize_channels_2(MDAIO_HEADER &H,MDAIO_HEADER &H_out,double *mean,double *stdev,int NN,FILE *input_file,FILE *output_file);

bool normalize_channels(const char *input_path,const char *output_path) {
	printf("normalize_channels %s %s...\n",input_path,output_path);
	FILE *input_file=fopen(input_path,"rb");
	if (!input_file) {
		printf("Unable to open file for reading: %s\n",input_path);
		return false;
	}
	FILE *output_file=fopen(output_path,"wb");
	if (!output_file) {
		fclose(input_file);
		printf("Unable to open file for writing: %s\n",output_path);
		return false;
	}
	MDAIO_HEADER H;
	mda_read_header(&H,input_file);
	int M=H.dims[0];
	int N=H.dims[1];

	int chunk_size=pow(2,15); // 2^15 = 32,768
	if (chunk_size>N) chunk_size=N;

	MDAIO_HEADER H_out=H;
	H_out.data_type=MDAIO_TYPE_FLOAT32;
	H_out.num_bytes_per_entry=4;
	mda_write_header(&H_out,output_file);

	QTime timer; timer.start();
	double *sum=(double *)malloc(sizeof(double)*M);
	double *sumsqr=(double *)malloc(sizeof(double)*M);
	double *mean=(double *)malloc(sizeof(double)*M);
	double *stdev=(double *)malloc(sizeof(double)*M);
	for (int m=0; m<M; m++) {
		sum[m]=0;
		sumsqr[m]=0;
	}
	for (int tt=0; tt<N; tt+=chunk_size) {
		if ((tt==0)||(timer.elapsed()>1000)) {
			float pct=tt*1.0/N;
			printf("Processing chunk %d/%d (%d%%)\n",(tt/chunk_size)+1,(N/chunk_size)+1,(int)(pct*100));
			timer.restart();
		}
		int tt2=tt+chunk_size;
		if (tt2>N) tt2=N;
		normalize_channels_1(H,sum,sumsqr,tt2-tt,input_file);
	}
	for (int m=0; m<M; m++) {
		mean[m]=sum[m]/N;
		//sum((a_i-mu)^2)=sumsqr - 2*(sum/N)*sum + N*(sum/N)^2 = sumsqr - 2*sum^2/N + sum^2/N = sumsqr-sum^2/N
		if (N>1) stdev[m]=sqrt((sumsqr[m]-sum[m]*sum[m]/N)/(N-1));
		else stdev[m]=1;
		if (stdev[m]==0) stdev[m]=1;
	}
	fseek(input_file,H.header_size,SEEK_SET);
	for (int tt=0; tt<N; tt+=chunk_size) {
		if ((tt==0)||(timer.elapsed()>1000)) {
			float pct=tt*1.0/N;
			printf("Scaling chunk %d/%d (%d%%)\n",(tt/chunk_size)+1,(N/chunk_size)+1,(int)(pct*100));
			timer.restart();
		}
		int tt2=tt+chunk_size;
		if (tt2>N) tt2=N;
		normalize_channels_2(H,H_out,mean,stdev,tt2-tt,input_file,output_file);
	}

	free(sum);
	free(sumsqr);
	free(mean);
	free(stdev);

	fclose(input_file);
	fclose(output_file);
	printf("[done].\n");

	return true;
}

void normalize_channels_1(MDAIO_HEADER &H,double *sum,double *sumsqr,int NN,FILE *input_file) {
	int M=H.dims[0];
	float *data=(float *)malloc(sizeof(float)*M*NN);
	mda_read_float32(data,&H,NN*M,input_file);
	int ii=0;
	for (int n=0; n<NN; n++) {
		for (int m=0; m<M; m++) {
			sum[m]+=data[ii+m];
			sumsqr[m]+=data[ii+m]*data[ii+m];
		}
		ii+=M;
	}
	free(data);
}

void normalize_channels_2(MDAIO_HEADER &H,MDAIO_HEADER &H_out,double *mean,double *stdev,int NN,FILE *input_file,FILE *output_file) {
	int M=H.dims[0];
	float *data=(float *)malloc(sizeof(float)*M*NN);
	mda_read_float32(data,&H,NN*M,input_file);
	int ii=0;
	for (int n=0; n<NN; n++) {
		for (int m=0; m<M; m++) {
			float val=data[ii+m];
			if (stdev[m]) {
				val=(val-mean[m])/stdev[m];
			}
			data[ii+m]=val;
		}
		ii+=M;
	}
	mda_write_float32(data,&H_out,NN*M,output_file);
	free(data);
}
