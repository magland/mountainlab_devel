#include "extract_clips.h"
#include "diskreadmda.h"
#include "mdaio.h"

int compute_K(DiskReadMda &X);

bool extract_clips(const char *input_path,const char *cluster_path,const char *output_path,const char *index_out_path,int clip_size) {
	DiskReadMda X(input_path);
	DiskReadMda C(cluster_path);

	if (X.totalSize()<=1) {
		printf("Problem reading input file: %s\n",input_path);
		return false;
	}
	if (C.totalSize()<=1) {
		printf("Problem reading input file: %s\n",cluster_path);
		return false;
	}

	int M=X.N1();
	int T=clip_size;
	int num_clips=C.N2();
	int K=compute_K(C);
	printf("K=%d\n",K);

	Mda index_out; index_out.allocate(1,K);

	MDAIO_HEADER H_out;
	H_out.data_type=MDAIO_TYPE_FLOAT32;
	H_out.num_bytes_per_entry=4;
	H_out.num_dims=3;
	H_out.dims[0]=M;
	H_out.dims[1]=T;
	H_out.dims[2]=num_clips;

	FILE *outf=fopen(output_path,"wb");
	if (!outf) {
		printf("Unable to open output file: %s\n",output_path);
		return false;
	}

	mda_write_header(&H_out,outf);

	float *buf=(float *)malloc(sizeof(float)*M*T);
	int jj=0;
	for (int k=1; k<=K; k++) {
		index_out.setValue(jj,0,k-1);
		for (int i=0; i<num_clips; i++) {
			int ii=0;
			int time0=(int)C.value(1,i);
			int k0=(int)C.value(2,i);
			if (k0==k) {
				for (int t=0; t<T; t++) {
					for (int m=0; m<M; m++) {
						buf[ii]=X.value(m,t+time0-T/2);
						ii++;
					}
				}
				mda_write_float32(buf,&H_out,M*T,outf);
				jj++;
			}
		}
	}
	free(buf);

	fclose(outf);

	if (!index_out.write(index_out_path)) {
		printf("Unable to write output file: %s\n",index_out_path);
		return false;
	}

	return true;
}

int compute_K(DiskReadMda &X) {
	int ret=0;
	for (int i=0; i<X.N2(); i++) {
		int k0=(int)X.value(2,i);
		if (k0>ret) ret=k0;
	}
	return ret;
}
