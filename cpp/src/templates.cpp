#include "templates.h"
#include "mdaio.h"
#include "mda.h"

bool compute_template(int M,int T,float *ret,const char *input_path,const Mda &CC,int k);
void extract_clip(int M,int T,float *ret,MDAIO_HEADER *H,FILE *input_file,int time0);

bool templates(const char *input_path,const char *cluster_path,const char *output_path,int clip_size) {

    printf("templates %s %s %s %d...\n",input_path,cluster_path,output_path,clip_size);


    Mda CC; CC.read(cluster_path);
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

    float *buf=(float *)malloc(sizeof(float)*M*clip_size);
    for (int k=1; k<=K; k++) {
        printf("Computing template %d/%d...\n",k,K);
        if (!compute_template(M,clip_size,buf,input_path,CC,k)) {
            fclose(output_file);
            return false;
        }
        mda_write_float32(buf,&H_out,M*clip_size,output_file);
    }
    free(buf);

    fclose(output_file);

    return true;
}

bool compute_template(int M,int T,float *ret,const char *input_path,const Mda &CC,int k) {

    FILE *input_file=fopen(input_path,"rb");
    if (!input_file) {
        printf("Unable to open file for reading: %s\n",input_path);
        return false;
    }

    MDAIO_HEADER H;
    mda_read_header(&H,input_file);

    double *sum=(double *)malloc(sizeof(double)*M*T);
    float *buf=(float *)malloc(sizeof(float)*M*T);
    for (int ii=0; ii<M*T; ii++) {
        sum[ii]=0;
    }
    int count=0;
    int NT=CC.N2();
    for (int ii=0; ii<NT; ii++) {
        int k0=(int)CC.value(2,ii);
        if (k0==k) {
            int time0=(int)CC.value(1,ii);
            extract_clip(M,T,buf,&H,input_file,time0);
            for (int ii=0; ii<M*T; ii++) {
                sum[ii]+=buf[ii];
            }
            count++;
        }
    }

    for (int ii=0; ii<M*T; ii++) {
        ret[ii]=sum[ii]/count;
    }

    free(sum);
    free(buf);

    fclose(input_file);

    return true;
}
void extract_clip(int M,int T,float *ret,MDAIO_HEADER *H,FILE *input_file,int time0) {
    for (int ii=0; ii<M*T; ii++) ret[ii]=0;
    int t1=time0-T/2;
    int t2=time0-T/2+T-1;
    int offset=0;
    if (t1<0) {
        offset=-t1; t1=0;
    }
    if (t2>=H->dims[1]) t2=H->dims[1]-1;
    fseek(input_file,H->header_size+H->num_bytes_per_entry*t1*M,SEEK_SET);
    mda_read_float32(&ret[offset*M],H,M*(t2-t1+1),input_file);
}
