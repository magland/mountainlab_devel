#include "remove_artifacts.h"

#include "diskreadmda.h"
#include "mdaio.h"

#include <QTime>

bool remove_artifacts(const char *input_path, const char *output_path, float threshold, bool normalize, int exclude_interval)
{
    DiskReadMda X(input_path);

    if (X.totalSize()<=1) {
        printf("Problem reading input file: %s\n",input_path);
        return false;
    }

    int M=X.N1();
    int N=X.N2();

    int *to_include=(int *)malloc(sizeof(int)*N);
    for (int i=0; i<N; i++) to_include[i]=1;
    float stdevs[M];
    for (int m=0; m<M; m++) stdevs[m]=1;

	QTime timer;
	timer.start();
	if (normalize) {
		printf("Normalizing...\n");
		int stride=1+N/10000;
		float sum[M]; float sumsqr[M]; int ct[M];
		for (int m=0; m<M; m++) {
			sum[m]=0; sumsqr[m]=0; ct[m]=0;
		}
		for (int i=0; i<N; i+=stride) {
			if (timer.elapsed()>1000) {
				printf("Part 1 - %d%%\n",(int)(i*1.0/N*100));
				timer.restart();
			}
			for (int m=0; m<M; m++) {
				float val=X.value(m,i);
				sum[m]+=val;
				sumsqr[m]+=val*val;
				ct[m]++;
			}
		}
		for (int m=0; m<M; m++) {
			if (ct[m]>=2) {
				stdevs[m]=sqrt((sumsqr[m]-sum[m]*sum[m]/ct[m])/(ct[m]-1));
			}
		}
	}

	printf("Identifying timepoints to remove...\n");
	timer.start();
	int i=0;
	while (i<N) {
		if ((i%1000==0)&&(timer.elapsed()>1000)) {
			printf("Part 2 - %d%%\n",(int)(i*1.0/N*100));
			timer.restart();
		}
		float maxval=0;
		for (int m=0; m<M; m++) {
			float val=fabs(X.value(m,i)/stdevs[m]);
			if (val>maxval) {
				maxval=val;
			}
		}
		if (maxval>threshold) {
			for (int dd=-exclude_interval/2; dd<-exclude_interval/2+exclude_interval; dd++) {
				int j=i+dd;
				if ((j>=0)&&(j<N)) {
					to_include[j]=0;
				}
			}
			i+=exclude_interval/6; //skip a bit to save some time
			i++;
		}
		else {
			i++;
		}
	}

    printf("Removing and writing output...\n");
    MDAIO_HEADER H_out;
    H_out.data_type=MDAIO_TYPE_FLOAT32;
    H_out.num_bytes_per_entry=4;
    H_out.num_dims=2;
    H_out.dims[0]=M;
    H_out.dims[1]=N;

    FILE *outf=fopen(output_path,"wb");
    if (!outf) {
        printf("Unable to open output file: %s\n",output_path);
        return false;
    }

    mda_write_header(&H_out,outf);

    int num_removed=0;
    float *buf=(float *)malloc(sizeof(float)*M);
    timer.start();
    for (int n=0; n<N; n++) {
		if ((n%1000==0)&&(timer.elapsed()>1000)) {
            printf("Part 3 - %d%%\n",(int)(n*1.0/N*100));
            timer.restart();
        }
        if (to_include[n]) {
            for (int m=0; m<M; m++) {
                buf[m]=X.value(m,n);
            }
        }
        else {
            num_removed++;
            for (int m=0; m<M; m++) {
                buf[m]=0;
            }
        }
        mda_write_float32(buf,&H_out,M,outf);
    }
    free(buf);

    fclose(outf);

    free(to_include);

    printf("%d timepoints removed (%g%%)",num_removed,num_removed*1.0/N*100);

    return true;
}
