#include "whiten.h"
#include "diskreadmda.h"
#include "get_principal_components.h"
#include "mdaio.h"

#include <QTime>
#include "Eigen/Core"
#include "Eigen/SVD"

void subtract_components(int M,int N,int ncomp,float *data,float *components);
void compute_svd(Mda &U,Mda &singvals,Mda &V,Mda &X);

bool whiten(const char *input_path,const char *output_path) {

	using namespace Eigen;

	DiskReadMda X;
	X.setPath(input_path);

	int M=X.N1();
	int N=X.N2();

	printf("Computing mean...\n");
	double mean[M];
	for (int m=0; m<M; m++) mean[m]=0;
	for (int n=0; n<N; n++) {
		for (int m=0; m<M; m++) {
			mean[m]+=X.value(m,n);
		}
	}
	for (int m=0; m<M; m++) {
		mean[m]/=N;
	}

	int N0=N;
	if (N0>1e6) N0=1e6;
	int stride=N/N0;

	printf("Setting up matrix Y...\n");
	MatrixXf Y(M,N0);
	int ii=0;
	for (int n0=0; n0<N0; n0++) {
		int n=n0*stride;
		for (int m=0; m<M; m++) {
			Y.data()[ii]=X.value(m,n)-mean[m];
			ii++;
		}
	}

	printf("Computing SVD...\n");
	JacobiSVD<MatrixXf> svd(Y, ComputeThinU | ComputeThinV);
	MatrixXf singvals0=svd.singularValues();
	MatrixXf U0=svd.matrixU();
	MatrixXf V0=svd.matrixV();

	MatrixXf D(M,M);
	for (int m2=0; m2<M; m2++) {
		for (int m1=0; m1<M; m1++) {
			if (m1==m2)
				D.data()[m1+M*m2]=1.0/singvals0.data()[m1];
			else
				D.data()[m1+M*m2]=0;
		}
	}

	printf("Setting matrix A...\n");
	MatrixXf A(M,M);
	//Somehow the following functionality is not supported by the version of libstdc++ that comes with MATLAB
	//A=U0*D*Transpose<MatrixXf>(U0);
	//so we need to do it by hand!
	for (int m2=0; m2<M; m2++) {
		for (int m1=0; m1<M; m1++) {
			float tmp=0;
			for (int i=0; i<M; i++) {
				tmp+=U0.data()[m1+M*i]*D.data()[i+M*i]*U0.data()[m2+M*i];
			}
			A.data()[m1+M*m2]=tmp;
		}
	}

	printf("Opening output file...\n");
	FILE *output_file=fopen(output_path,"wb");
	if (!output_file) {
		printf("Unable to open file for writing: %s\n",output_path);
		return false;
	}

	MDAIO_HEADER H_out;
	H_out.data_type=MDAIO_TYPE_FLOAT32;
	H_out.num_bytes_per_entry=4;
	H_out.num_dims=2;
	H_out.dims[0]=M;
	H_out.dims[1]=N;
	mda_write_header(&H_out,output_file);

	printf("Whitening and writing output...\n");
	for (int n=0; n<N; n++) {
		float input[M];
		for (int m=0; m<M; m++) {
			input[m]=X.value(m,n)-mean[m];
		}
		float output[M];
		for (int m=0; m<M; m++) output[m]=0;
		int kk=0;
		for (int mc=0; mc<M; mc++) {
			float coeff=input[mc];
			for (int mr=0; mr<M; mr++) {
				output[mr]+=A.data()[kk]*coeff;
				kk++;
			}
		}
		mda_write_float32(output,&H_out,M,output_file);
	}

	printf("Closing file.\n");
	fclose(output_file);

	return true;
}

void compute_svd(Mda &U,Mda &singvals,Mda &V,Mda &X) {
	using namespace Eigen;
	MatrixXf X0(X.N1(),X.N2());
	int ii=0;
	for (int c=0; c<X.N2(); c++) {
		for (int r=0; r<X.N1(); r++) {
			X0.data()[ii]=X.value(r,c);
			ii++;
		}
	}

	JacobiSVD<MatrixXf> svd(X0, ComputeThinU | ComputeThinV);
	MatrixXf singvals0=svd.singularValues();
	MatrixXf U0=svd.matrixU();
	MatrixXf V0=svd.matrixV();

	singvals.allocate(singvals0.rows(),singvals0.cols());
	ii=0;
	for (int c=0; c<singvals.N2(); c++) {
		for (int r=0; r<singvals.N1(); r++) {
			singvals.setValue(singvals0.data()[ii],r,c);
			ii++;
		}
	}

	U.allocate(U0.rows(),U0.cols());
	ii=0;
	for (int c=0; c<U.N2(); c++) {
		for (int r=0; r<U.N1(); r++) {
			U.setValue(U0.data()[ii],r,c);
			ii++;
		}
	}

	V.allocate(V0.rows(),V0.cols());
	ii=0;
	for (int c=0; c<V.N2(); c++) {
		for (int r=0; r<V.N1(); r++) {
			V.setValue(V0.data()[ii],r,c);
			ii++;
		}
	}
}

bool whiten_old(const char *input_path,const char *output_path,int ncomp) {

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
