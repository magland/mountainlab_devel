#include "cross_correlograms.h"
#include "diskreadmda.h"
#include <QList>
#include <stdio.h>

typedef QList<int> IntList;
bool cross_correlograms(const char *clusters_path,const char *output_path,int max_dt) {
	QList<int> times,labels;

    printf("Setting up times and labels...\n");
	DiskReadMda C; C.setPath(clusters_path);
	int L=C.N2();
	int K=1;
	for (int ii=0; ii<L; ii++) {
		int time0=(int)C.value(1,ii);
		int label0=(int)C.value(2,ii);
		times << time0;
		labels << label0;
		if (label0>K) K=label0;
	}

    printf("Initializing output...\n");
	QList<IntList> out;
	IntList empty_list;
	for (int k1=1; k1<=K; k1++) {
		for (int k2=1; k2<=K; k2++) {
			out << empty_list;
		}
	}

    printf("Setting time differences...\n");
	int i1=0;
	for (int i2=0; i2<L; i2++) {
        while ((i1+1<L)&&(times[i1]<times[i2]-max_dt)) i1++;
		int k2=labels[i2];
		int t2=times[i2];
        if (k2>=1) {
            for (int jj=i1; jj<i2; jj++) {
                int k1=labels[jj];
                int t1=times[jj];
                if (k1>=1) {
                    out[(k1-1)+K*(k2-1)] << t2-t1;
                    out[(k2-1)+K*(k1-1)] << t1-t2;
                }
            }
        }
	}

    printf("Counting...\n");
	int ct=0;
	for (int k1=1; k1<=K; k1++) {
		for (int k2=1; k2<=K; k2++) {
			ct+=out[(k1-1)+K*(k2-1)].count();
		}
	}

    printf("Creating mda...\n");
	Mda ret; ret.allocate(3,ct);

	ct=0;
	for (int k1=1; k1<=K; k1++) {
		for (int k2=1; k2<=K; k2++) {
			IntList *tmp=&out[(k1-1)+K*(k2-1)];
			for (int jj=0; jj<tmp->count(); jj++) {
				ret.setValue(k1,0,ct);
				ret.setValue(k2,1,ct);
				ret.setValue((*tmp)[jj],2,ct);
				ct++;
			}
		}
	}

    printf("Writing...\n");
	ret.write(output_path);

	return true;

}


