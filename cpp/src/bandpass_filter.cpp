#include "bandpass_filter.h"

#include "mdaio.h"
#include "fftw3.h"
#include <QTime>
#include <math.h>
#include "omp.h"

bool do_fft_1d_r2c(int N,float *out,float *in);
bool do_ifft_1d_c2r(int N,float *out,float *in);
void multiply_complex_by_real_kernel(int N,float *Y,float *kernel);
void define_kernel(int N,float *kernel,double samplefreq,double freq_min,double freq_max);
void bandpass_filter_2(MDAIO_HEADER &H,FILE *input_file,MDAIO_HEADER &H_out,FILE *output_file,long timepoint1,long timepoint2,double samplefreq,double freq_min,double freq_max,long overlap,double outlier_threshold);

bool bandpass_filter(const char *input_path,const char *output_path,double samplefreq,double freq_min,double freq_max,double outlier_threshold) {
	printf("bandpass filter %s %s %g %g %g...\n",input_path,output_path,samplefreq,freq_min,freq_max);
	FILE *input_file=fopen(input_path,"rb");
	if (!input_file) {
		printf("Unable to open file for reading: %s\n",input_path);
		return false;
	}
	printf("%s...",output_path);
	FILE *output_file=fopen(output_path,"wb");
	if (!output_file) {
		fclose(input_file);
		printf("Unable to open file for writing: %s\n",output_path);
		return false;
	}
	MDAIO_HEADER H;
	mda_read_header(&H,input_file);
	//int M=H.dims[0];
	int N=H.dims[1];

    long chunk_size=pow(2,12); // 2^15 = 32,768
	if (chunk_size>N) chunk_size=N;
	long overlap=chunk_size/10;
	if (overlap<2000) overlap=2000;

	MDAIO_HEADER H_out=H;
	H_out.data_type=MDAIO_TYPE_FLOAT32;
	H_out.num_bytes_per_entry=4;
	mda_write_header(&H_out,output_file);

	QTime timer; timer.start();
	for (long tt=0; tt<N; tt+=chunk_size) {
		if ((tt==0)||(timer.elapsed()>1000)) {
			float pct=tt*1.0/N;
			printf("Processing chunk %ld/%ld (%d%%)\n",(tt/chunk_size)+1,(N/chunk_size)+1,(int)(pct*100));
			timer.restart();
		}
		long tt2=tt+chunk_size;
		if (tt2>N) tt2=N;
		bandpass_filter_2(H,input_file,H_out,output_file,tt,tt2,samplefreq,freq_min,freq_max,overlap,outlier_threshold);
	}

	fclose(input_file);
	fclose(output_file);
	printf("[done].\n");

	return true;
}

void bandpass_filter_2(MDAIO_HEADER &H,FILE *input_file,MDAIO_HEADER &H_out,FILE *output_file,long timepoint1,long timepoint2,double samplefreq,double freq_min,double freq_max,long overlap,double outlier_threshold) {
	int M=H.dims[0];
	long N=timepoint2-timepoint1;

	long overlap1=overlap;
	if (timepoint1-overlap1<0) overlap1=timepoint1;
	long overlap2=overlap;
	if (timepoint2+overlap2>H.dims[1]) overlap2=H.dims[1]-timepoint2;
	long NN=N+overlap1+overlap2;

	float *kernel0=(float *)malloc(sizeof(float)*NN);
	define_kernel(NN,kernel0,samplefreq,freq_min,freq_max);

	float *X=(float *)malloc(sizeof(float)*M*NN);
	fseek(input_file,H.header_size+H.num_bytes_per_entry*M*(timepoint1-overlap1),SEEK_SET);
	mda_read_float32(X,&H,M*NN,input_file);

    int num_threads=M;
    omp_set_num_threads(num_threads);
    #pragma omp parallel
	{

		double *X0=(double *)fftw_malloc(sizeof(double)*NN);
		fftw_complex *Y0=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*NN);
		fftw_plan p_fft;
		#pragma omp critical
		p_fft=fftw_plan_dft_r2c_1d(NN,X0,Y0,FFTW_ESTIMATE);
		fftw_plan p_ifft;
		#pragma omp critical
		p_ifft=fftw_plan_dft_c2r_1d(NN,Y0,X0,FFTW_ESTIMATE);

        //if (omp_get_thread_num()==0) printf("Using %d threads.\n",omp_get_num_threads());
        #pragma omp for
		for (int m=0; m<M; m++) {
			int ii=m;
			for (long n=0; n<NN; n++) {
				X0[n]=X[ii];
				ii+=M;
			}
			if (outlier_threshold>0) {
				for (long n=0; n<NN; n++) {
					if (X0[n]<-outlier_threshold) X0[n]=-outlier_threshold;
					if (X0[n]>outlier_threshold) X0[n]=outlier_threshold;
				}
			}
			fftw_execute(p_fft);
			for (long n=0; n<NN; n++) {
				Y0[n][0]*=kernel0[n];
				Y0[n][1]*=kernel0[n];
			}
			fftw_execute(p_ifft);
			ii=m;
			for (long n=0; n<NN; n++) {
				X[ii]=X0[n]/NN; ii+=M;
			}
		}

		fftw_free(X0);
		fftw_free(Y0);
	}

	mda_write_float32(&X[M*overlap1],&H_out,M*N,output_file);

	free(X);
	free(kernel0);
}

bool do_ifft_1d_c2r(int N,float *out,float *in) {
	/*
	if (num_threads>1) {
		fftw_init_threads();
		fftw_plan_with_nthreads(num_threads);
	}
	*/

	fftw_complex *in2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex *out2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	for (int ii=0; ii<N; ii++) {
		in2[ii][0]=in[ii*2];
		in2[ii][1]=in[ii*2+1];
	}
	fftw_plan p;
	#pragma omp critical
	p=fftw_plan_dft_1d(N,in2,out2,FFTW_BACKWARD,FFTW_ESTIMATE);

	fftw_execute(p);
	for (int ii=0; ii<N; ii++) {
		//out[ii*2]=out2[ii][0];
		//out[ii*2+1]=out2[ii][1];
		out[ii]=out2[ii][0];
	}
	fftw_free(in2);
	fftw_free(out2);

	/*
	if (num_threads>1) {
		fftw_cleanup_threads();
	}
	*/

	#pragma omp critical
	fftw_destroy_plan(p);

	return true;
}


bool do_fft_1d_r2c(int N,float *out,float *in) {
	/*
	if (num_threads>1) {
		fftw_init_threads();
		fftw_plan_with_nthreads(num_threads);
	}
	*/

	fftw_complex *in2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex *out2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N);
	for (int ii=0; ii<N; ii++) {
		//in2[ii][0]=in[ii*2];
		//in2[ii][1]=in[ii*2+1];
		in2[ii][0]=in[ii];
		in2[ii][1]=0;
	}
	fftw_plan p;
	#pragma omp critical
	p=fftw_plan_dft_1d(N,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);

	fftw_execute(p);
	for (int ii=0; ii<N; ii++) {
		out[ii*2]=out2[ii][0];
		out[ii*2+1]=out2[ii][1];
	}
	fftw_free(in2);
	fftw_free(out2);

	/*
	if (num_threads>1) {
		fftw_cleanup_threads();
	}
	*/

	#pragma omp critical
	fftw_destroy_plan(p);

	return true;
}
void multiply_complex_by_real_kernel(int N,float *Y,float *kernel) {
	for (int i=0; i<N; i++) {
		Y[i*2]*=kernel[i];
		Y[i*2+1]*=kernel[i];
	}
}

void define_kernel(int N,float *kernel,double samplefreq,double freq_min,double freq_max) {
	//Based on ahb's MATLAB code
	float T=N/samplefreq; //total time
	//frequency grid
	float df=1/T;
	float *fgrid=(float *)malloc(sizeof(float)*N);
	for (int i=0; i<N; i++) {
		if (i<=(N+1)/2) fgrid[i]=df*i;
		else fgrid[i]=df*(i-N);
	}

    float fwidlo=100; // roll-off width (Hz). Sets ringing timescale << 10 ms
    float fwidhi=1000; // roll-off width (Hz). Sets ringing timescale << 1 ms

	for (int i=0; i<N; i++) {
		float absf=fabs(fgrid[i]);
		float val=1;
		val*=(1+tanh((absf-freq_min)/fwidlo))/2;
		val*=(1-tanh((absf-freq_max)/fwidhi))/2;
		kernel[i]=val;
	}

	free(fgrid);
}
