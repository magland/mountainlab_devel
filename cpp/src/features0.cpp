#include "features.h"
#include "mdaio.h"
#include "pcasolver.h"
#include "get_principal_components.h"
#include "diskreadmda.h"
#include <QDebug>

int features_2(int ch,DiskReadMda &X,DiskReadMda &D,DiskReadMda &A,MDAIO_HEADER *H_out,FILE *output_file,int num_features,int clip_size);

bool features(const char *input_path,const char *detect_path,const char *adjacency_path,const char *output_path,int num_features,int clip_size) {
    printf("features %s %s %s %d %d...\n",input_path,detect_path,output_path,num_features,clip_size);
    FILE *output_file=fopen(output_path,"wb");
    if (!output_file) {
        printf("Unable to open file for writing: %s\n",output_path);
        return false;
    }

    DiskReadMda X; X.setPath(input_path);
    DiskReadMda D; D.setPath(detect_path);
    DiskReadMda A; A.setPath(adjacency_path);

    int M=X.N1();

    // channel, timepoint, feature1, feature2, ..., feature6
    MDAIO_HEADER H_out;
    H_out.data_type=MDAIO_TYPE_FLOAT32;
    H_out.num_bytes_per_entry=4;
    H_out.num_dims=2;
    H_out.dims[0]=num_features+2;
    H_out.dims[1]=2; //to be replaced later
    mda_write_header(&H_out,output_file);

    int total_num_events=0;
	for (int ch=1; ch<=M; ch++) {
		printf("ch=%d... ",ch);
		int num_events=features_2(ch,X,D,A,&H_out,output_file,num_features,clip_size);
		printf("%d events...\n",num_events);
		total_num_events+=num_events;
	}

    H_out.dims[1]=total_num_events;
    fseek(output_file,0,SEEK_SET);
    mda_write_header(&H_out,output_file);

    fclose(output_file);

    return true;
}

int features_2(int ch,DiskReadMda &X,DiskReadMda &D,DiskReadMda &A,MDAIO_HEADER *H_out,FILE *output_file,int num_features,int clip_size) {
    int M=X.N1();
    int N=X.N2();
    int NT=D.N2();

    int M0=0;
    for (int m=0; m<M; m++) {
        if (A.value(m,ch-1)) {
            M0++;
        }
    }
    int channels[M0];
    int mm=0;
    for (int m=0; m<M; m++) {
        if (A.value(m,ch-1)) {
            channels[mm]=m;
            mm++;
        }
    }

    int num_events=0;
    for (int i=0; i<NT; i++) {
        if (D.value(0,i)==ch) {
            int time0=(int)D.value(1,i);
            if ((time0-clip_size/2>=0)&&(time0-clip_size/2+clip_size<=N)) {
                num_events++;
            }
        }
    }
    if (num_events==0) return 0;
    float *Y=(float *)malloc(sizeof(float)*M0*clip_size*num_events);
    int *times=(int *)malloc(sizeof(int)*num_events);
    int ie=0;
	for (int i=0; i<NT; i++) {
		if (D.value(0,i)==ch) {
			int time0=(int)D.value(1,i);
			if ((time0-clip_size/2>=0)&&(time0-clip_size/2+clip_size<=N)) {
				int aa=M0*clip_size*ie;
				for (int t=0; t<clip_size; t++) {
					for (int im=0; im<M0; im++) {
						Y[aa]=X.value(channels[im],t+time0-clip_size/2);
						aa++;
					}
				}
				times[ie]=time0;
				ie++;
			}
		}
	}

    /*
    float *components=(float *)malloc(sizeof(float)*M0*clip_size*num_features);
    get_principal_components(M0*clip_size,num_events,num_features,components,Y);
    for (int ie=0; ie<num_events; ie++) {
        float features[num_features];
        int time0=times[ie];
        for (int cc=0; cc<num_features; cc++) {
            float *component=&components[M0*clip_size*cc];
            float ip=0;
            for (int t=0; t<clip_size; t++) {
                for (int im=0; im<M; im++) {
                    ip+=component[im+M0*t]*Y[im+M0*t+M0*clip_size*ie];
                }
            }
            features[cc]=ip;
        }
        float cc_tt[2]; cc_tt[0]=ch; cc_tt[1]=time0;
        mda_write_float32(cc_tt,H_out,2,output_file);
        mda_write_float32(features,H_out,num_features,output_file);
    }
    */
    float *all_features=(float *)malloc(sizeof(float)*num_features*N);
	get_pca_features(M0*clip_size,num_events,num_features,all_features,Y);
    for (int ie=0; ie<num_events; ie++) {
        int time0=times[ie];
        float cc_tt[2]; cc_tt[0]=ch; cc_tt[1]=time0;
        mda_write_float32(cc_tt,H_out,2,output_file);
        mda_write_float32(&all_features[num_features*ie],H_out,num_features,output_file);
    }

    free(all_features);

    free(Y);
    //free(components);

    return num_events;
}

