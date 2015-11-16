#include "extract.h"
#include "mdaio.h"

bool extract(const char *input_path,const char *output_path,int num_channels,int M,int *channels,long t1,long t2,float threshold) {
	printf("extract data %s %s...\n",input_path,output_path);
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
	MDAIO_HEADER H_out;
	H_out.num_dims=2;
	H_out.dims[0]=M;
	H_out.dims[1]=t2-t1;
	H_out.data_type=MDAIO_TYPE_FLOAT32;
	H_out.num_bytes_per_entry=4;
	mda_write_header(&H_out,output_file);

	short int *buf=(short int *)malloc(sizeof(short int)*num_channels);
	float *buf2=(float *)malloc(sizeof(float)*M);
	for (long t=t1; t<t2; t++) {
		int num0=fread(buf,sizeof(short int),num_channels,input_file);
		if (num0!=num_channels) {
			printf("Problem reading raw input file (t=%ld/%ld).\n",t,t2);
			fclose(input_file);
			fclose(output_file);
			return false;
		}
		for (int m=0; m<M; m++) {
			float val=buf[channels[m]];
			if (threshold>0) {
				if (val<-threshold) val=-threshold;
				if (val>threshold) val=threshold;
			}
			buf2[m]=val;
		}
		mda_write_float32(buf2,&H_out,M,output_file);
	}
	free(buf);

	fclose(input_file);
	fclose(output_file);

	return true;
}
