#include "detect.h"
#include <stdio.h>
#include "mdaio.h"
#include <math.h>
#include <QTime>

int detect_2(MDAIO_HEADER &H,FILE *input_file,MDAIO_HEADER &H_out,FILE *output_file,long timepoint1,long timepoint2,long overlap,int inner_window_width,int outer_window_width,double threshold,bool normalize,bool individual_channels);

bool detect(const char *input_path,const char *output_path,int inner_window_width,int outer_window_width,double threshold,bool normalize,bool individual_channels,int clip_size) {
	printf("detect %s %s %d %d %g...\n",input_path,output_path,inner_window_width,outer_window_width,threshold);
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

	long chunk_size=pow(2,15); // 2^15 = 32,768
	if (chunk_size>N) chunk_size=N;
	//long overlap=chunk_size/10;
	long overlap=chunk_size;
	if (overlap<2000) overlap=2000;

	MDAIO_HEADER H_out;
	H_out.data_type=MDAIO_TYPE_INT32;
	H_out.num_bytes_per_entry=4;
	H_out.num_dims=2;
	H_out.dims[0]=2;
	H_out.dims[1]=2; //to be determined later
	mda_write_header(&H_out,output_file);

	QTime timer; timer.start();
	int total_num_events=0;
	for (long tt=0; tt<N; tt+=chunk_size) {
		if ((tt==0)||(timer.elapsed()>1000)) {
			float pct=tt*1.0/N;
			printf("Processing chunk %ld/%ld (%d%%)\n",(tt/chunk_size)+1,(N/chunk_size)+1,(int)(pct*100));
			timer.restart();
		}
		long tt1=tt;
		long tt2=tt+chunk_size;
		if (tt1<clip_size) tt1=clip_size;
		if (tt2>N-clip_size) tt2=N-clip_size;
		total_num_events+=detect_2(H,input_file,H_out,output_file,tt1,tt2,overlap,inner_window_width,outer_window_width,threshold,normalize,individual_channels);
	}
	printf("total_num_events=%d\n",total_num_events);
	fseek(output_file,0,SEEK_SET);
	H_out.dims[1]=total_num_events;
	mda_write_header(&H_out,output_file); //now that H_out.dims[1] has been set, rewrite the header

	fclose(input_file);
	fclose(output_file);
	printf("[done].\n");

	return true;
}

int detect_2(MDAIO_HEADER &H,FILE *input_file,MDAIO_HEADER &H_out,FILE *output_file,long timepoint1,long timepoint2,long overlap,int inner_window_width,int outer_window_width,double threshold,bool normalize,bool individual_channels) {
	int M=H.dims[0];
	long N=timepoint2-timepoint1;

	long overlap1=overlap;
	if (timepoint1-overlap1<0) overlap1=timepoint1;
	long overlap2=overlap;
	if (timepoint2+overlap2>H.dims[1]) overlap2=H.dims[1]-timepoint2;
	long NN=N+overlap1+overlap2;

	float *X=(float *)malloc(sizeof(float)*M*NN);
	fseek(input_file,H.header_size+H.num_bytes_per_entry*M*(timepoint1-overlap1),SEEK_SET);
	mda_read_float32(X,&H,M*NN,input_file);

	int *tmp_output=(int *)malloc(sizeof(int)*M*NN);

	if (individual_channels) {
		#pragma omp parallel
		{

			double *X0=(double *)malloc(sizeof(double)*NN);
			int *detected_inds1=(int *)malloc(sizeof(int)*NN);
			int *use_it=(int *)malloc(sizeof(int)*NN);

			//if (omp_get_thread_num()==0) printf("Using %d threads.\n",omp_get_num_threads());
			#pragma omp for
			for (int m=0; m<M; m++) {
				int num_detected_inds=0;
				for (int n=0; n<NN; n++) {
					tmp_output[m+n*M]=-1;
				}
				int ii=m;
				for (long n=0; n<NN; n++) {
					X0[n]=X[ii];
					ii+=M;
				}
				double sliding_sumsqr=0;
				int sliding_count=0;
                if (normalize) {
                    for (int n=0; (n<=outer_window_width/2)&&(n<NN); n++) {
                        sliding_sumsqr+=X0[n]*X0[n];
                        sliding_count++;
                    }
                }
				for (long n=0; n<NN; n++) {
                    if (normalize) {
                        if (n+outer_window_width/2<NN) {
                            float val0=X0[n+outer_window_width/2];
                            sliding_sumsqr+=val0*val0;
                            sliding_count++;
                        }
                        if (n-outer_window_width/2-1>=0) {
                            float val0=X0[n-outer_window_width/2-1];
                            sliding_sumsqr-=val0*val0;
                            sliding_count--;
                        }
                    }
					if (n>inner_window_width) {
						float val00=X0[n];
                        double stdev=1;
                        if (normalize) stdev=sqrt(sliding_sumsqr/(sliding_count-1));
						if (stdev>0) {
							float absval00=val00; if (val00<0) absval00=-val00;
							if (absval00/stdev>=threshold) {
								detected_inds1[num_detected_inds]=n;
								num_detected_inds++;
							}
						}
					}
				}
				for (long i=0; i<num_detected_inds; i++) {
					use_it[i]=1;
					float val=X0[detected_inds1[i]];
					int k=i-1;
					while ((k>=0)&&(detected_inds1[k]>=detected_inds1[i]-inner_window_width)) {
						if (fabs(val)>fabs(X0[detected_inds1[k]])) use_it[k]=0;
						if (fabs(val)<fabs(X0[detected_inds1[k]])) use_it[i]=0;
						k--;
					}
				}
				int jjj=0;
				for (long i=0; i<num_detected_inds; i++) {
					if (use_it[i]) {
						tmp_output[m+jjj*M]=detected_inds1[i];
						jjj++;
					}
				}
			}

			free(detected_inds1);
			free(use_it);
			free(X0);
		}
	}
	else {
		#pragma omp parallel
		{

			double *X0=(double *)malloc(sizeof(double)*NN);
			int *detected_inds1=(int *)malloc(sizeof(int)*NN);
			int *use_it=(int *)malloc(sizeof(int)*NN);

			//if (omp_get_thread_num()==0) printf("Using %d threads.\n",omp_get_num_threads());
			//#pragma omp for
			{
				int num_detected_inds=0;
				for (int n=0; n<NN; n++) {
					for (int m=0; m<M; m++) {
						tmp_output[m+n*M]=-1;
					}
				}
				int ii=0;
				for (long n=0; n<NN; n++) {
					float best_val=0;
					for (int m=0; m<M; m++) {
						if (fabs(X[ii])>=fabs(best_val)) {
							best_val=X[ii];
						}
						ii++;
					}
					X0[n]=best_val;
				}
				double sliding_sumsqr=0;
				int sliding_count=0;
                if (normalize) {
                    for (int n=0; (n<=outer_window_width/2)&&(n<NN); n++) {
                        sliding_sumsqr+=X0[n]*X0[n];
                        sliding_count++;
                    }
                }
				for (long n=0; n<NN; n++) {
                    if (normalize) {
                        if (n+outer_window_width/2<NN) {
                            float val0=X0[n+outer_window_width/2];
                            sliding_sumsqr+=val0*val0;
                            sliding_count++;
                        }
                        if (n-outer_window_width/2-1>=0) {
                            float val0=X0[n-outer_window_width/2-1];
                            sliding_sumsqr-=val0*val0;
                            sliding_count--;
                        }
                    }
					if (n>inner_window_width) {
						float val00=X0[n];
                        double stdev=1;
                        if (normalize) stdev=sqrt(sliding_sumsqr/(sliding_count-1));
						if (stdev>0) {
							float absval00=val00; if (val00<0) absval00=-val00;
							if (absval00/stdev>=threshold) {
								detected_inds1[num_detected_inds]=n;
								num_detected_inds++;
							}
						}
					}
				}
				for (long i=0; i<num_detected_inds; i++) {
					use_it[i]=1;
					float val=X0[detected_inds1[i]];
					int k=i-1;
					while ((k>=0)&&(detected_inds1[k]>=detected_inds1[i]-inner_window_width)) {
						if (fabs(val)>fabs(X0[detected_inds1[k]])) use_it[k]=0;
						if (fabs(val)<fabs(X0[detected_inds1[k]])) use_it[i]=0;
						k--;
					}
				}
				int jjj=0;
				for (long i=0; i<num_detected_inds; i++) {
					if (use_it[i]) {
						int m0=0; //always associate with first channel!
						tmp_output[m0+jjj*M]=detected_inds1[i];
						jjj++;
					}
				}
			}

			free(detected_inds1);
			free(use_it);
			free(X0);
		}
	}

	int total_num_events=0;
	for (int m=0; m<M; m++) {
		int i=0;
		while ((i<NN)&&(tmp_output[m+i*M]>=0)) {
			int n=tmp_output[m+i*M];
			if ((n>=overlap1)&&(n<overlap1+N)) { //it's important to get this right!
				total_num_events++;
			}
			i++;
		}
	}
	int *Y=(int *)malloc(sizeof(int)*2*total_num_events);
	int ie=0;
	for (int m=0; m<M; m++) {
		int i=0;
		while ((i<NN)&&(tmp_output[m+i*M]>=0)) {
			int n=tmp_output[m+i*M];
			if ((n>=overlap1)&&(n<overlap1+N)) { //it's important to get this right!
				Y[0+2*ie]=m +1;
				Y[1+2*ie]=timepoint1+n-overlap1 +1;
				ie++;
			}
			i++;
		}
	}

	int ct=mda_write_int32(Y,&H_out,2*total_num_events,output_file);
	if (ct!=2*total_num_events) {
		printf("WARNING: unexpected problem writing output! %d<>%d",ct,2*total_num_events);
	}

	free(tmp_output);
	free(X);
	free(Y);

	return total_num_events;
}
