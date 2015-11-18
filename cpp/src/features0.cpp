#include "features.h"
#include "mdaio.h"
#include "pcasolver.h"
#include "get_principal_components.h"

bool get_header(MDAIO_HEADER *H,const char *path);
int features_2(int ch,const char *input_path,const char *detect_path,MDAIO_HEADER *H_out,FILE *output_file,int num_features,int clip_size);
float compute_inner_product(int N,float *X,float *Y);

bool features(const char *input_path,const char *detect_path,const char *output_path,int num_features,int clip_size) {
    printf("features %s %s %s %d %d...\n",input_path,detect_path,output_path,num_features,clip_size);
    FILE *output_file=fopen(output_path,"wb");
    if (!output_file) {
        printf("Unable to open file for writing: %s\n",output_path);
        return false;
    }

    MDAIO_HEADER H_input;
    if (!get_header(&H_input,input_path)) {
        printf("Error reading %s\n",input_path);
        fclose(output_file); return false;
    }

    MDAIO_HEADER H_detect;
    if (!get_header(&H_detect,detect_path)) {
        printf("Error reading %s\n",detect_path);
        fclose(output_file); return false;
    }

    int M=H_input.dims[0];

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
        int num_events=features_2(ch,input_path,detect_path,&H_out,output_file,num_features,clip_size);
        printf("%d events...\n",num_events);
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

int features_2(int ch,const char *input_path,const char *detect_path,MDAIO_HEADER *H_out,FILE *output_file,int num_features,int clip_size) {


    FILE *input_file=fopen(input_path,"rb");
    if (!input_file) {
        printf("Unable to open file for reading: %s\n",input_path);
        return -1;
    }
    FILE *detect_file=fopen(detect_path,"rb");
    if (!detect_file) {
        fclose(input_file);
        printf("Unable to open file for reading: %s\n",detect_path);
        return -1;
    }

    MDAIO_HEADER H_input;
    mda_read_header(&H_input,input_file);
    MDAIO_HEADER H_detect;
    mda_read_header(&H_detect,detect_file);

    int M=H_input.dims[0];
    int NT=H_detect.dims[1];

    float *X=0; // IN THE FUTURE WE CAN USE A SUBSET OF THE DATA TO COMPUTE THE PCA
    float *components=0;
    float *buf=(float *)malloc(sizeof(float)*M*clip_size);
    int N=0;

    for (int pass=1; pass<=3; pass++) {
        fseek(detect_file,H_detect.header_size,SEEK_SET);
        int ii=0;
        for (int i=0; i<NT; i++) {
            float tmp[2];
            mda_read_float32(tmp,&H_detect,2,detect_file);
            if (tmp[0]==ch) {
                int time0=(int)tmp[1];
                if ((time0-clip_size/2>=0)&&(time0-clip_size/2+clip_size<=H_input.dims[1])) {
                    if (pass==2) {
                        fseek(input_file,H_input.header_size+(M*time0-clip_size/2)*H_input.num_bytes_per_entry,SEEK_SET);
                        mda_read_float32(&X[M*clip_size*ii],&H_input,M*clip_size,input_file);
                    }
                    if (pass==3) {
                        fseek(input_file,H_input.header_size+(M*time0-clip_size/2)*H_input.num_bytes_per_entry,SEEK_SET);
                        mda_read_float32(buf,&H_input,M*clip_size,input_file);
                        float features[num_features];
                        for (int cc=0; cc<num_features; cc++) {
                            features[cc]=compute_inner_product(M*clip_size,buf,&components[M*clip_size*cc]);
                        }
                        float cc_tt[2]; cc_tt[0]=ch; cc_tt[1]=time0;
                        mda_write_float32(cc_tt,H_out,2,output_file);
                        mda_write_float32(features,H_out,num_features,output_file);
                    }
                    ii++;
                }
            }
        }
        if (pass==1) {
            N=ii;
            if (N>0) X=(float *)malloc(sizeof(float)*M*clip_size*N);
        }
        if (pass==2) {
            if (N>0) {
                // IN THE FUTURE WE CAN USE A SUBSET OF THE DATA TO COMPUTE THE PCA
                components=(float *)malloc(sizeof(float)*M*clip_size*num_features);
                get_principal_components(M*clip_size,N,num_features,components,X);
            }
        }
    }

    free(buf);
    if (components) free(components);
    if (X) free(X);

    fclose(input_file);
    fclose(detect_file);

    return N;
}

bool get_header(MDAIO_HEADER *H,const char *path) {
    FILE *file=fopen(path,"rb");
    if (!file) {
        printf("Unable to open file for reading num channels: %s\n",path);
        return false;
    }
    mda_read_header(H,file);
    fclose(file);
    return true;
}

float compute_inner_product(int N,float *X,float *Y) {
    float ret=0;
    for (int i=0; i<N; i++) {
        ret+=X[i]*Y[i];
    }
    return ret;
}
