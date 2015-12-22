#include "fit.h"
#include <stdio.h>
#include "mdaio.h"
#include <math.h>
#include <QTime>
#include <QList>
#include "mda.h"

int fit_2(MDAIO_HEADER &H,FILE *input_file,MDAIO_HEADER &H_out,FILE *output_file,long timepoint1,long timepoint2,long overlap,int inner_window_width,int outer_window_width,double threshold);

struct Event {
    int channel;
    int timepoint;
    int label;
};

bool fit(const char *input_path,const char *templates_path,const char *cluster_in_path,const char *cluster_out_path) {
    printf("fit %s %s %s %s..\n",input_path,templates_path,cluster_in_path,cluster_out_path);

    FILE *input_file=fopen(input_path,"rb");
    if (!input_file) {
        printf("Unable to open file for reading: %s\n",input_path);
        return false;
    }

    Mda templates; templates.read(templates_path);
    Mda cluster_in; cluster_in.read(cluster_in_path);
    QList<Event> events_in;
    for (int i=0; i<cluster_in.N2(); i++) {
        Event EE;
        EE.channel=(int)cluster_in.value(0,i);
        EE.timepoint=(int)cluster_in.value(1,i);
        EE.label=(int)cluster_in.value(2,i);
        events_in << EE;
    }

    /*
    long chunk_size=pow(2,15); // 2^15 = 32,768
    if (chunk_size>N) chunk_size=N;
    long overlap=chunk_size;
    if (overlap<2000) overlap=2000;

    QTime timer; timer.start();
    int total_num_events=0;
    for (long tt=0; tt<N; tt+=chunk_size) {
        if ((tt==0)||(timer.elapsed()>1000)) {
            float pct=tt*1.0/N;
            printf("Processing chunk %ld/%ld (%d%%)\n",(tt/chunk_size)+1,(N/chunk_size)+1,(int)(pct*100));
            timer.restart();
        }
        long tt2=tt+chunk_size;
        if (tt2>N) tt2=N;

        int timepoint1=tt;
        int timepoint2=tt2;

        long overlap1=overlap;
        if (timepoint1-overlap1<0) overlap1=timepoint1;
        long overlap2=overlap;
        if (timepoint2+overlap2>H.dims[1]) overlap2=H.dims[1]-timepoint2;
        long NN=N+overlap1+overlap2;

        float *X=(float *)malloc(sizeof(float)*M*NN);
        fseek(input_file,H.header_size+H.num_bytes_per_entry*M*(timepoint1-overlap1),SEEK_SET);
        mda_read_float32(X,&H,M*NN,input_file);

        fit_2(M,NN,X,templates,cluster_in,cluster_out);

        free(X);
    }
    printf("total_num_events=%d\n",total_num_events);
    fseek(output_file,0,SEEK_SET);
    H_out.dims[1]=total_num_events;
    mda_write_header(&H_out,output_file); //now that H_out.dims[1] has been set, rewrite the header

    fclose(input_file);
    fclose(output_file);
    printf("[done].\n");
    */

    return true;
}

