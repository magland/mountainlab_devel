#include "features.h"
#include "mdaio.h"
#include "pcasolver.h"
#include "get_principal_components.h"
#include "diskreadmda.h"
#include <QDebug>
#include "mda.h"

Mda features_2(int ch,DiskReadMda &X,DiskReadMda &D,const Mda &A,int num_features,int clip_size);

bool features(const char *input_path,const char *detect_path,const char *adjacency_path,const char *output_path,int num_features,int clip_size) {
    printf("features %s %s %s %d %d...\n",input_path,detect_path,output_path,num_features,clip_size);

	Mda A; A.read(adjacency_path);

	int M=0;
	{
		DiskReadMda X; X.setPath(input_path);
		M=X.N1();
	}

	QList<Mda> all_features;
	for (int ch=1; ch<=M; ch++) {
		Mda X;
		all_features << X;
	}

	#pragma omp parallel
	{
		#pragma omp for
		for (int ch=1; ch<=M; ch++) {
			DiskReadMda X; X.setPath(input_path); //important to do these here? for multi-threaded safety?
			DiskReadMda D; D.setPath(detect_path);
			printf("ch=%d... ",ch);
			Mda features=features_2(ch,X,D,A,num_features,clip_size);
			int num_events=features.N2();
			if (num_events==1) num_events=0; //handle the terrible case because mda can't be of size zero -- :(
			printf("%d events...\n",num_events);
			all_features[ch-1]=features;
		}
	}

	int total_num_events=0;
	for (int ch=0; ch<M; ch++) {
		int num0=all_features[ch].N2();
		if (num0==1) num0=0; //handle the terrible case because mda can't be of size zero -- :(
		total_num_events+=all_features[ch].N2();
	}

	Mda output; output.allocate(num_features+2,total_num_events);
	int ie=0;
	for (int ch=0; ch<M; ch++) {
		Mda X=all_features[ch];
		if (X.totalSize()>1) { //handle the terrible case because mda can't be of size zero -- :(
			for (int ii=0; ii<X.N2(); ii++) {
				for (int jj=0; jj<X.N1(); jj++) {
					output.setValue(X.value(jj,ii),jj,ie);
				}
				ie++;
			}
		}
	}

	output.write(output_path);

    return true;
}

Mda features_2(int ch,DiskReadMda &X,DiskReadMda &D,const Mda &A,int num_features,int clip_size) {
    int M=X.N1();
    int N=X.N2();
    int NT=D.N2();

	Mda empty_mda; empty_mda.allocate(1,1); //problem, we are not allowed to have an mda of size 0 -- :(

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
            int time0=(int)D.value(1,i)-1; //convert to zero-based indexing
            if ((time0-clip_size/2>=0)&&(time0-clip_size/2+clip_size<=N)) {
                num_events++;
            }
        }
    }
	if (num_events==0) return empty_mda;
    float *Y=(float *)malloc(sizeof(float)*M0*clip_size*num_events);
    int *times=(int *)malloc(sizeof(int)*num_events);
    int ie=0;
	for (int i=0; i<NT; i++) {
		if (D.value(0,i)==ch) {
            int time0=(int)D.value(1,i)-1; //convert to zero-based indexing
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

    float *all_features=(float *)malloc(sizeof(float)*num_features*N);
	get_pca_features(M0*clip_size,num_events,num_features,all_features,Y);

	Mda ret;
	ret.allocate(num_features+2,num_events);
    for (int ie=0; ie<num_events; ie++) {
		int time0=times[ie];
		ret.setValue(ch,0,ie);
        ret.setValue(time0+1,1,ie); //convert back to 1-based indexing for the output
		for (int ff=0; ff<num_features; ff++) {
			ret.setValue(all_features[num_features*ie+ff],ff+2,ie);
		}
    }

    free(all_features);
    free(Y);
    //free(components);

	return ret;
}

