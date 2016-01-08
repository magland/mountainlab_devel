#include "templates.h"
#include "mdaio.h"
#include "diskreadmda.h"

bool compute_template(int M,int T,float *ret,DiskReadMda &X,DiskReadMda &CC,int k);
void extract_clip(int M,int T,float *ret,MDAIO_HEADER *H,FILE *input_file,int time0);

bool templates(const char *input_path,const char *cluster_path,const char *output_path,int clip_size) {

    printf("templates %s %s %s %d...\n",input_path,cluster_path,output_path,clip_size);

    DiskReadMda CC; CC.setPath(cluster_path);
    int NT=CC.N2();

    int K=0,M=0;
    for (int ii=0; ii<NT; ii++) {
        int ch=(int)CC.value(0,ii);
        int k=(int)CC.value(2,ii);
        if (ch>M) M=ch;
        if (k>K) K=k;
    }

    if (K==0) {
        printf("Problem: K is zero\n");
        return false;
    }
    if (M==0) {
        printf("Problem: M is zero\n");
        return false;
    }

    FILE *output_file=fopen(output_path,"wb");
    if (!output_file) {
        printf("Unable to open file for writing: %s\n",output_path);
        return false;
    }

    MDAIO_HEADER H_out;
    H_out.data_type=MDAIO_TYPE_FLOAT32;
    H_out.num_bytes_per_entry=4;
    H_out.num_dims=3;
    H_out.dims[0]=M;
    H_out.dims[1]=clip_size;
    H_out.dims[2]=K;
    mda_write_header(&H_out,output_file);

    DiskReadMda X; X.setPath(input_path);

    float *buf=(float *)malloc(sizeof(float)*M*clip_size);
	printf("Computing %d templates ...\n",K);
    for (int k=1; k<=K; k++) {
        if (!compute_template(M,clip_size,buf,X,CC,k)) {
            fclose(output_file);
            return false;
        }
        mda_write_float32(buf,&H_out,M*clip_size,output_file);
    }
    free(buf);

    fclose(output_file);

    return true;
}

bool compute_template(int M,int T,float *ret,DiskReadMda &X,DiskReadMda &CC,int k) {
    double *sum=(double *)malloc(sizeof(double)*M*T);
    for (int ii=0; ii<M*T; ii++) {
        sum[ii]=0;
    }
    int count=0;
    int NT=CC.N2();
    for (int ii=0; ii<NT; ii++) {
        int k0=(int)CC.value(2,ii);
        if (k0==k) {
            int time0=(int)CC.value(1,ii)-1;
            if ((time0-T/2>=0)&&(time0-T/2+T<=X.N2())) {
                int jj=0;
                for (int t=time0-T/2; t<time0-T/2+T; t++) {
                    for (int m=0; m<M; m++) {
                        sum[jj]+=X.value(m,t);
                        jj++;
                    }
                }
                count++;
            }
        }
    }

    if (count>0) {
        for (int ii=0; ii<M*T; ii++) {
            ret[ii]=sum[ii]/count;
        }
    }

    free(sum);

    return true;
}
