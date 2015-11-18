#include "cluster.h"
#include "mdaio.h"
#include "isosplit.h"

int get_num_channels(const char *path);
int cluster_2(int ch,const char *input_path,MDAIO_HEADER *H_out,FILE *output_file,int *num_clusters);
int do_cluster(int M,int N,float *X,int *labels);

bool cluster(const char *input_path,const char *output_path) {
    printf("cluster %s %s...\n",input_path,output_path);
    FILE *output_file=fopen(output_path,"wb");
    if (!output_file) {
        printf("Unable to open file for writing: %s\n",output_path);
        return false;
    }

    int M=get_num_channels(input_path);

    // channel, timepoint, label
    MDAIO_HEADER H_out;
    H_out.data_type=MDAIO_TYPE_FLOAT32;
    H_out.num_bytes_per_entry=4;
    H_out.num_dims=2;
    H_out.dims[0]=3;
    H_out.dims[1]=2; //to be replaced later
    mda_write_header(&H_out,output_file);

    int total_num_events=0;
    for (int ch=1; ch<=M; ch++) {
        printf("ch=%d... ",ch);
        int num_clusters;
        int num_events=cluster_2(ch,input_path,&H_out,output_file,&num_clusters);
        printf("%d events, %d clusters...\n",num_events,num_clusters);
        if (num_events<0) {
            fclose(output_file);
            return false;
        }
        total_num_events+=num_events;
    }

    H_out.dims[1]=total_num_events;
    fseek(output_file,0,SEEK_SET);
    mda_write_header(&H_out,output_file);

    fclose(output_file);

    return true;
}

int cluster_2(int ch,const char *input_path,MDAIO_HEADER *H_out,FILE *output_file,int *num_clusters) {
    FILE *input_file=fopen(input_path,"rb");
    if (!input_file) {
        printf("Unable to open file for reading: %s\n",input_path);
        return -1;
    }

    MDAIO_HEADER H_input;
    mda_read_header(&H_input,input_file);

    int num_features=H_input.dims[0]-2;
    int NT=H_input.dims[1];

    float *X=0;
    int *times=0;
    int N=0;
    for (int pass=1; pass<=2; pass++) {
        fseek(input_file,H_input.header_size,SEEK_SET);
        int ii=0;
        for (int i=0; i<NT; i++) {
            float tmp[2+num_features];
            mda_read_float32(tmp,&H_input,2+num_features,input_file);
            if (tmp[0]==ch) {
                int time0=(int)tmp[1];
                if (pass==2) {
                    for (int j=0; j<num_features; j++) {
                        X[num_features*ii+j]=tmp[2+j];
                    }
                    times[ii]=time0;
                }
                ii++;
            }
        }
        if (pass==1) {
            N=ii;
            if (N>0) {
                X=(float *)malloc(sizeof(float)*num_features*N);
                times=(int *)malloc(sizeof(int)*N);
            }
        }
    }
    if (N>0) {
        int *labels=(int *)malloc(sizeof(int)*N);
        *num_clusters=do_cluster(num_features,N,X,labels);
        for (int n=0; n<N; n++) {
            float buf[3];
            buf[0]=ch;
            buf[1]=times[n];
            buf[2]=labels[n];
            mda_write_float32(buf,H_out,3,output_file);
        }
        free(labels);
    }

    if (times) free(times);
    if (X) free(X);

    fclose(input_file);

    return N;
}

int get_num_channels(const char *path) {
    int ret=0;

    FILE *file=fopen(path,"rb");
    if (!file) {
        printf("Unable to open file for getting num channels: %s\n",path);
        return 0;
    }
    MDAIO_HEADER H;
    mda_read_header(&H,file);

    int num_features=H.dims[0]-2;
    int NT=H.dims[1];
    float buf[num_features+2];
    for (int ii=0; ii<NT; ii++) {
        mda_read_float32(buf,&H,num_features+2,file);
        int ch=(int)buf[0];
        if (ch>ret) ret=ch;
    }

    fclose(file);
    return ret;
}

int do_cluster(int M,int N,float *X,int *labels) {
    Mda A; A.allocate(M,N);
    int ii=0;
    for (int n=0; n<N; n++) {
        for (int m=0; m<M; m++) {
            A.setValue(X[ii],m,n);
            ii++;
        }
    }
    QVector<int> labels0=isosplit(A);
    int ret=0;
    for (int n=0; n<N; n++) {
        labels[n]=labels0.value(n);
        if (labels[n]>ret) ret=labels[n];
    }
    return ret;
}
