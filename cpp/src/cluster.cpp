#include "cluster.h"
#include "mdaio.h"
#include "isosplit.h"
#include "diskreadmda.h"

int get_num_channels(const char *path);
int cluster_2(int ch,DiskReadMda &F,MDAIO_HEADER *H_out,FILE *output_file,int *num_clusters,int label_offset);
int do_cluster(int M,int N,float *X,int *labels);

bool cluster(const char *input_path,const char *output_path) {
    printf("cluster %s %s...\n",input_path,output_path);
    FILE *output_file=fopen(output_path,"wb");
    if (!output_file) {
        printf("Unable to open file for writing: %s\n",output_path);
        return false;
    }

    int M=get_num_channels(input_path);
    DiskReadMda F; F.setPath(input_path);

    // channel, timepoint, label
    MDAIO_HEADER H_out;
    H_out.data_type=MDAIO_TYPE_FLOAT32;
    H_out.num_bytes_per_entry=4;
    H_out.num_dims=2;
    H_out.dims[0]=3;
    H_out.dims[1]=2; //to be replaced later
    mda_write_header(&H_out,output_file);

    int total_num_events=0;
    int label_offset=0;
    for (int ch=1; ch<=M; ch++) {
        printf("ch=%d... ",ch);
        int num_clusters;
        int num_events=cluster_2(ch,F,&H_out,output_file,&num_clusters,label_offset);
        label_offset+=num_clusters;
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

int cluster_2(int ch,DiskReadMda &F,MDAIO_HEADER *H_out,FILE *output_file,int *num_clusters,int label_offset) {
    int num_features=F.N1()-2;
    int NT=F.N2();

    int num_events=0;
    for (int i=0; i<NT; i++) {
        if (F.value(0,i)==ch) {
            num_events++;
        }
    }
    if (num_events==0) {
        *num_clusters=0;
        return 0;
    }

    float *X=(float *)malloc(sizeof(float)*num_features*num_events);
    int *times=(int *)malloc(sizeof(int)*num_events);

    int ii=0;
    int ie=0;
    for (int i=0; i<NT; i++) {
        if (F.value(0,i)==ch) {
            for (int f=0; f<num_features; f++) {
                X[ii]=F.value(f+2,i);
                ii++;
            }
            times[ie]=(int)F.value(1,i);
            ie++;
        }
    }

    int *labels=(int *)malloc(sizeof(int)*num_events);
    *num_clusters=do_cluster(num_features,num_events,X,labels);
    for (int ie=0; ie<num_events; ie++) {
        float buf[3];
        buf[0]=ch;
        buf[1]=times[ie];
        buf[2]=labels[ie]+label_offset;
        mda_write_float32(buf,H_out,3,output_file);
    }

    //debug
    {
        Mda debug; debug.allocate(num_features,num_events);
        int jjj=0;
        for (int ie=0; ie<num_events; ie++) {
            for (int f=0; f<num_features; f++) {
                debug.setValue(X[jjj],f,ie);
                jjj++;
            }
        }
        debug.write(QString("debug-ch-%1.mda").arg(ch));
    }
    {
        Mda debug; debug.allocate(1,num_events);
        for (int ie=0; ie<num_events; ie++) {
            debug.setValue(labels[ie],0,ie);
        }
        debug.write(QString("debug-labels-ch-%1.mda").arg(ch));
    }
    //////////////////////////////////////////////

    free(labels);

    free(X);
    free(times);

    return num_events;
}

int get_num_channels(const char *path) {
    DiskReadMda X; X.setPath(path);
    int ret=0;
    int NT=X.N2();
    for (int ii=0; ii<NT; ii++) {
        int ch=X.value(0,ii);
        if (ch>ret) ret=ch;
    }
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
