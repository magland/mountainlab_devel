#include "extract_channels.h"

#include "diskreadmda.h"
#include "mdaio.h"

bool extract_channels(const char *input_path, const char *output_path, const QList<int> &channels)
{
    DiskReadMda X(input_path);

    if (X.totalSize()<=1) {
        printf("Problem reading input file: %s\n",input_path);
        return false;
    }

    //int M1=X.N1();
    int N=X.N2();

    int M2=channels.count();

    MDAIO_HEADER H_out;
    H_out.data_type=MDAIO_TYPE_FLOAT32;
    H_out.num_bytes_per_entry=4;
    H_out.num_dims=2;
    H_out.dims[0]=M2;
    H_out.dims[1]=N;

    FILE *outf=fopen(output_path,"wb");
    if (!outf) {
        printf("Unable to open output file: %s\n",output_path);
        return false;
    }

    mda_write_header(&H_out,outf);

    float *buf=(float *)malloc(sizeof(float)*M2);
    for (int n=0; n<N; n++) {
        for (int m2=0; m2<M2; m2++) {
            buf[m2]=X.value(channels[m2]-1,n); //convert to zero-based indexing
        }
        mda_write_float32(buf,&H_out,M2,outf);
    }
    free(buf);

    fclose(outf);

    return true;
}
