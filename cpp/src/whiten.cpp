#include "whiten.h"
#include "mdaio.h"
#include "get_principal_components.h"

#include <QTime>

void subtract_components(int M,int N,int ncomp,float *data,float *components);

bool whiten(const char *input_path,const char *output_path,int ncomp) {

	printf("Opening input file %s...\n",input_path);
	FILE *input_file=fopen(input_path,"rb");
	if (!input_file) {
		printf("Unable to open file for reading: %s\n",input_path);
		return false;
	}
	printf("Opening output file %s...\n",output_path);
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

	MDAIO_HEADER H_out=H;
	H_out.data_type=MDAIO_TYPE_FLOAT32;
	H_out.num_bytes_per_entry=4;
	mda_write_header(&H_out,output_file);

	int NN=10000; if (NN>N) NN=N;
	float *data=(float *)malloc(sizeof(float)*M*NN);
	int stride=N/NN;

	printf("Reading data...\n");
	for (int i=0; i<NN; i++) {
		mda_read_float32(&data[M*i],&H,M,input_file);
		fseek(input_file,M*(stride-1)*H.num_bytes_per_entry,SEEK_CUR);
	}

	float *components=(float *)malloc(sizeof(float)*M*ncomp);
	printf("Computing principal components...\n");
	//get_principal_components_2(M,NN,ncomp,components,data);
	if (ncomp>0) get_pca_components(M,NN,ncomp,components,data);
	fseek(input_file,H.header_size,SEEK_SET);
	printf("Subtracting %d components...\n",ncomp);
	QTime timer; timer.start();
	int chunk_size=10000;
	float *data0=(float *)malloc(sizeof(float)*M*chunk_size);
	for (int i=0; i<N; i+=chunk_size) {
		if ((i==0)||(i==N-1)||(timer.elapsed()>1000)) {
			float pct=(i+1)*1.0/N;
			printf("%d%%\n",(int)(pct*100));
			timer.restart();
		}
		int i2=i+chunk_size; if (i2>N) i2=N;
		mda_read_float32(data0,&H,M*(i2-i),input_file);
		if (ncomp>0) subtract_components(M,i2-i,ncomp,data0,components);
		mda_write_float32(data0,&H_out,M*(i2-i),output_file);
	}


    free(components);
	free(data);
	free(data0);

	fclose(input_file);
	fclose(output_file);

	printf("Done with whiten.\n");

	return true;
}

void subtract_components(int M,int N,int ncomp,float *data,float *components) {
	for (int n=0; n<N; n++) {
		float *data1=&data[M*n];
		for (int cc=0; cc<ncomp; cc++) {
			float *C=&components[cc*M];
			float ip=0;
			for (int m=0; m<M; m++) {
				ip+=C[m]*data1[m];
			}
			for (int m=0; m<M; m++) {
				data1[m]-=ip*C[m];
			}
		}
	}
}
