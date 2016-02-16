#include "split_firings.h"
#include "diskreadmda.h"
#include "mda.h"
#include "get_principal_components.h"
#include "isosplit.h"
#include "omp.h"
#include <stdio.h>

QList<int> get_times(DiskReadMda &C,int k);
int get_K(DiskReadMda &C);
void extract_clips(Mda &clips,DiskReadMda &X,QList<int> &times,int clip_size);
void compute_features(Mda &features,Mda &clips,int num_features);
int get_max_k(const QVector<int> &labels);

bool split_firings(const char *input_path,const char *firings_path,const char *output_path,int num_features,int clip_size,float ks_threshold,int k_init) {

	DiskReadMda C; C.setPath(firings_path);
    Mda C2;
    C2.allocate(C.N1(),C.N2());
    for (int i=0; i<C.N2(); i++) {
        C2.setValue(C.value(0,i),0,i);
        C2.setValue(C.value(1,i),1,i);
    }

    int K=get_K(C);

    int kk=1;
	//#pragma omp parallel for (is this a good idea?)
    for (int k=1; k<=K; k++) {
		DiskReadMda X; X.setPath(input_path);  //needed here for thread safety?
		DiskReadMda C; C.setPath(firings_path);
		printf("k=%d/%d... ",k,K);
        QList<int> times=get_times(C,k);
        Mda clips;
		printf("extract clips... ");
		extract_clips(clips,X,times,clip_size);
        Mda features;
		printf("compute features... ");
		compute_features(features,clips,num_features);
		printf("isosplit... ");
		QVector<int> labels0=isosplit(features,ks_threshold,k_init);
		printf("setting...\n");
        int K0=get_max_k(labels0);
		if (K0>1) printf("::: split into %d clusters\n",K0);
		else printf("\n");
		int jjj=0;
		for (int bb=0; bb<C.N2(); bb++) {
			if (C.value(2,bb)==k) {
			   C2.setValue(kk+labels0[jjj]-1,2,bb);
			   jjj++;
			}
		}
		kk+=K0;
    }

	C2.write(output_path);

	return true;
}

int get_max_k(const QVector<int> &labels) {
	int ret=0;
	for (int i=0; i<labels.count(); i++) {
		if (labels[i]>ret) ret=labels[i];
	}
	return ret;
}

QList<int> get_times(DiskReadMda &C,int k) {
	QList<int> times;
	for (int i=0; i<C.N2(); i++) {
		if (C.value(2,i)==k) {
            times << (int)C.value(1,i)-1; //convert to zero-based indexing
		}
	}
	return times;
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
	clips.allocate(M,clip_size,times.count());
	int tt1=-clip_size/2;
	for (int i=0; i<times.count(); i++) {
		int t0=times[i];
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

	get_pca_features(M*T,num,num_features,features.dataPtr(),clips.dataPtr());
}
