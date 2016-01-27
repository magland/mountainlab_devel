#include "extract_clips.h"
#include "diskreadmda.h"
#include "mdaio.h"

int compute_K(DiskReadMda &X);

bool extract_clips(const char *input_path,const char *detect_path,const char *output_clips_path,int clip_size) {
	DiskReadMda X(input_path);
    DiskReadMda D(detect_path);

	if (X.totalSize()<=1) {
		printf("Problem reading input file: %s\n",input_path);
		return false;
	}
    if (D.totalSize()<=1) {
        printf("Problem reading input file: %s\n",detect_path);
		return false;
	}

    int M=X.N1();
	int T=clip_size;
    int num_clips=D.N2();

	MDAIO_HEADER H_out;
	H_out.data_type=MDAIO_TYPE_FLOAT32;
	H_out.num_bytes_per_entry=4;
	H_out.num_dims=3;
	H_out.dims[0]=M;
	H_out.dims[1]=T;
	H_out.dims[2]=num_clips;

    FILE *outf=fopen(output_clips_path,"wb");
	if (!outf) {
        printf("Unable to open output file: %s\n",output_clips_path);
		return false;
	}

	mda_write_header(&H_out,outf);

	float *buf=(float *)malloc(sizeof(float)*M*T);
    for (int ii=0; ii<num_clips; ii++) {
        int time0=(int)D.value(1,ii)-1; //convert to zero-based indexing
        int jj=0;
        for (int t=0; t<T; t++) {
            for (int m=0; m<M; m++) {
                buf[jj]=X.value(m,t+time0-T/2);
                jj++;
            }
        }
        mda_write_float32(buf,&H_out,M*T,outf);
	}
	free(buf);

	fclose(outf);

	return true;
}
