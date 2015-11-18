#include "features.h"
//#include "mdaio.h"

int get_num_channels(const char *path);
/*
bool detect(const char *input_path,const char *detect_path,const char *output_path,int num_features,int clip_size) {
    printf("features %s %s %s %d %d...\n",input_path,detect_path,output_path,num_features,clip_size);
    printf("%s...",output_path);
    FILE *output_file=fopen(output_path,"wb");
    if (!output_file) {
        printf("Unable to open file for writing: %s\n",output_path);
        return false;
    }
    int M=get_num_channels(input_path);
    if (M==0) {
        printf("Problem reading number of channels.\n");
        return false;
    }

    // channel, timepoint, feature1, feature2, ..., feature6
    MDAIO_HEADER H_out;
    H_out.data_type=MDAIO_TYPE_FLOAT32;
    H_out.num_bytes_per_entry=4;
    H_out.num_dims=2;
    H_out.dims[0]=num_features+2;
    H_out.dims[1]=2;
    mda_write_header(&H_out,output_file);

    for (int m=0; m<M; m++) {
        if (!detect_2(input_path,detect_path,output_file,num_features,clip_size)) {
            fclose(output_file);
            return false;
        }
    }

    fclose(output_file);

    return true;
}

int get_num_channels(const char *path) {
//    FILE *file=fopen(output_path,"wb");
//    if (!file) {
//        printf("Unable to open file for reading num channels: %s\n",path);
//        return 0;
//    }
//    MDAIO_HEADER H;
//    mda_read_header(&H,file);
//    fclose(file);
//    return H.dims[0];
}
*/
