#include "split_clusters.h"
#include "diskreadmda.h"
#include "mda.h"
#include "get_principal_components.h"
#include "isosplit.h"

int get_K(DiskReadMda &C);
void extract_clips(Mda &clips,DiskReadMda &X,QList<int> &times,int clip_size);
void compute_features(Mda &features,Mda &clips,int num_features);

bool cluster(const char *input_path,const char *cluster_path,const char *output_path,int num_features,int clip_size) {
    DiskReadMda X; X.setPath(input_path);
    DiskReadMda C; C.setPath(cluster_path);
    Mda C2;
    C2.allocate(C.N1(),C.N2());
    for (int i=0; i<C.N2(); i++) {
        C2.setValue(C.value(0,i),0,i);
        C2.setValue(C.value(1,i),1,i);
    }

    int M=X.N1();
    int K=get_K(C);

    int kk=1;
    for (int k=1; k<=K; k++) {
        QList<int> times=get_times(C,k);
        Mda clips;
		extract_clips(clips,X,times,clip_size);
        Mda features;
		compute_features(features,clips,num_features);
        QVector<int> labels0=isosplit(features);
        int K0=get_max_k(labels0);
        for (aa=1; aa<=K0; aa++) {
            int jjj=0;
            for (int bb=0; bb<C.N2(); bb++) {
                if (C.value(2,bb)==k) {
                   C2.setValue(kk,2,jjj);
                   jjj++;
                }
            }
            kk++;
        }
    }
}

int get_K(DiskReadMda &C) {
    int ret=0;
    for (int i=0; i<C.N2(); i++) {
        int k0=C.value(2,i);
        if (k0>ret) ret=k0;
    }
    return ret;
}

void extract_clips(Mda &clips,DiskReadMda &X,QList<int> &times,int clip_size) {
	int M=X.N1();
	clips.allocate(M,clip_size,);
	int tt1=-clip_size/2;
	int tt2=tt1+clip_size-1;
	for (int i=0; i<times.count(); i++) {
		for (int t=0; t<clip_size; t++) {
			for (int m=0; m<M; m++) {
				clips.setValue(X.value(m,t0+t+tt1),m,t,i);
			}
		}
	}
}

void compute_features(Mda &features,Mda &clips,int num_features) {

	int M=clips.N1();
	int T=clips.N2();
	int num=clips.N3();
	features.allocate(num_features,num);
	get_principal_components(M*T,num,num_features,features.dataPtr(),clips.dataPtr());
}
