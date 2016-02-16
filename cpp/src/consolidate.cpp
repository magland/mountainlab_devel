#include "consolidate.h"
#include "diskreadmda.h"
#include "mdaio.h"
#include "mda.h"
#include <QSet>

//void get_load_channels(int K,int *load_channels,DiskReadMda &cluster);
Mda get_load_channels(int K,DiskReadMda &cluster);

struct SortRec {
	double val;
	int index;
};
bool caseInsensitiveLessThan(const SortRec &R1, const SortRec &R2)
{
  return R1.val<R2.val;
}
QList<int> get_sort_order(const QList<double> &values) {
	QList<SortRec> list;
	for (int i=0; i<values.count(); i++) {
		SortRec RR; RR.val=values[i]; RR.index=i;
		list << RR;
	}
	qSort(list.begin(),list.end(),caseInsensitiveLessThan);
	QList<int> ret;
	for (int i=0; i<list.count(); i++) {
		ret << list[i].index;
	}
	return ret;
}
QList<int> get_template_sort_order(const Mda &X) {
	int K=X.N3();
	QList<double> values;
	for (int k=0; k<K; k++) {
		double sum1=0;
		double sum2=0;
		for (int t=0; t<X.N2(); t++) {
			for (int m=0; m<X.N1(); m++) {
				double val0=X.value(m,t,k);
				sum1+=val0*val0*m;
				sum2+=val0*val0;
			}
		}
		if (sum2>0) sum1/=sum2;
		values << sum1;
	}
	QList<int> inds=get_sort_order(values);
	return inds;
}

typedef QSet<int> IntSet;

float compute_overlap(const QSet<int> &S1,const QSet<int> &S2) {
	if ((S1.count()==0)||(S2.count()==0)) return 0;
	QSet<int> S3=S1;
	S3.intersect(S2);
	int ct=S3.count();
	float ret1=ct*1.0/S1.count();
	float ret2=ct*1.0/S2.count();
	if (ret1>ret2) return ret1;
	else return ret2;
}

bool consolidate(const char *firings_path,const char *templates_path,const char *cluster_out_path,const char *templates_out_path,const char *load_channels_out_path,float coincidence_threshold) {
	DiskReadMda cluster; cluster.setPath(firings_path);
	Mda templates; templates.read(templates_path);

	int M=templates.N1();
	int T=templates.N2();
	int K=templates.N3();

	Mda load_channels=get_load_channels(K,cluster);


	bool to_use[K+1];
	for (int k=1; k<=K; k++) { //loop through each spike type
		float energies[M+1]; //compute the per-channel energies
		for (int m=1; m<=M; m++) {
			float sumsqr=0;
			for (int t=0; t<T; t++) {
				float val=templates.value(m-1,t,k-1);
				sumsqr+=val*val;
			}
			energies[m]=sumsqr;
		}
		int ch=load_channels.value(0,k-1);
		bool okay=true;
		for (int m=1; m<=M; m++) {
			if (m!=ch) {
				if (energies[ch]<=0.9*energies[m]) { //if the energy is higher on any of the other channels, remove it
					okay=false;
				}
			}
		}
		to_use[k]=okay;
	}

	{
		QList<IntSet> times_bin_20;
		IntSet empty;
		for (int k=0; k<=K; k++) {
			times_bin_20 << empty;
		}
		for (int i=0; i<cluster.N2(); i++) {
			int k0=(int)cluster.value(2,i);
			int t0=(int)cluster.value(1,i);
			if (k0<times_bin_20.count()) {
				times_bin_20[k0].insert(t0/20);
			}
		}
		for (int k1=1; k1<=K; k1++) {
			for (int k2=1; k2<=K; k2++) {
				if (k1!=k2) {
					if ((to_use[k1])&&(to_use[k2])) {
						float overlap_frac=compute_overlap(times_bin_20[k1],times_bin_20[k2]);
						if (overlap_frac>0.25) {
							printf("OVERLAP (%d%%) BETWEEN CLUSTERS %d and %d.\n",(int)(overlap_frac*100),k1,k2);
						}
						if (overlap_frac>coincidence_threshold) {
							if (times_bin_20[k1].count()<times_bin_20[k2].count()) {
								printf("SUFFICIENT OVERLAP (%d%%) BETWEEN CLUSTERS %d and %d, removing %d.\n",(int)(overlap_frac*100),k1,k2,k1);
								to_use[k1]=false;
							}
							else {
								printf("SUFFICIENT OVERLAP (%d%%) BETWEEN CLUSTERS %d and %d, removing %d.\n",(int)(overlap_frac*100),k1,k2,k2);
								to_use[k2]=false;
							}
						}
					}
				}
			}
		}
	}

	int num_to_use=0;
	for (int i=1; i<=K; i++) {
		if (to_use[i]) num_to_use++;
	}

	int K2=num_to_use;
	int num_events2=0;
	for (int i=0; i<cluster.N2(); i++) {
		int k=(int)cluster.value(2,i);
		if (to_use[k]) num_events2++;
	}

	int kk=1;
	int cluster_mapping1[K+1];
	for (int i=1; i<=K; i++) {
		if (to_use[i]) {
			cluster_mapping1[i]=kk;
			kk++;
		}
		else {
			cluster_mapping1[i]=0;
		}
	}

	Mda templates1;
	templates1.allocate(M,T,K2);
	for (int k=1; k<=K; k++) {
		if (cluster_mapping1[k]) {
			for (int t=0; t<T; t++) {
				for (int m=0; m<M; m++) {
					templates1.setValue(templates.value(m,t,k-1),m,t,cluster_mapping1[k]-1);
				}
			}
		}
	}
	QList<int> order=get_template_sort_order(templates1);

	int cluster_mapping2[K+1];
	for (int a=0; a<=K; a++) cluster_mapping2[a]=0;
	for (int k=1; k<=K2; k++) {
		int k2=order[k-1]+1;
		for (int a=1; a<=K; a++) {
			if (cluster_mapping1[a]==k2) {
				cluster_mapping2[a]=k;
			}
		}
	}
	for (int k=1; k<=K; k++) {
		if (cluster_mapping2[k]) {
			printf("MAPPING %d->%d\n",k,cluster_mapping2[k]);
		}
	}

	Mda load_channels2; load_channels2.allocate(1,K2);
	Mda templates2; templates2.allocate(M,T,K2);
	for (int k=1; k<=K; k++) {
		if (cluster_mapping2[k]) {
			load_channels2.setValue(load_channels.value(0,k-1),0,cluster_mapping2[k]-1);
			for (int t=0; t<T; t++) {
				for (int m=0; m<M; m++) {
					templates2.setValue(templates.value(m,t,k-1),m,t,cluster_mapping2[k]-1);
				}
			}
		}
	}

	MDAIO_HEADER H_cluster;
	H_cluster.data_type=MDAIO_TYPE_FLOAT32;
	H_cluster.num_bytes_per_entry=4;
	H_cluster.num_dims=2;
	H_cluster.dims[0]=3;
	H_cluster.dims[1]=num_events2;
	{
		FILE *outf=fopen(cluster_out_path,"wb");
		if (!outf) {
			printf("Unable to open file for writing: %s\n",cluster_out_path);
			return false;
		}
		mda_write_header(&H_cluster,outf);
		for (int i=0; i<cluster.N2(); i++) {
			int k0=(int)cluster.value(2,i);
			if (cluster_mapping2[k0]) {
				int k=cluster_mapping2[k0];
				float ch_tt_k[3];
				ch_tt_k[0]=cluster.value(0,i);
				ch_tt_k[1]=cluster.value(1,i);
				ch_tt_k[2]=k;
				mda_write_float32(ch_tt_k,&H_cluster,3,outf);
			}
		}
		fclose(outf);
	}


	templates2.write(templates_out_path);
	load_channels2.write(load_channels_out_path);

	printf("Using %d clusters\n",num_to_use);

	return true;
}

Mda get_load_channels(int K,DiskReadMda &cluster) {
	Mda load_channels;
	load_channels.allocate(1,K);
	for (int k=0; k<K; k++) load_channels.setValue(0,0,k);
	for (int i=0; i<cluster.N2(); i++) {
		int ch=(int)cluster.value(0,i);
		int k=(int)cluster.value(2,i);
		if ((1<=k)&&(k<=K)) load_channels.setValue(ch,0,k-1);
	}
	return load_channels;
}

/*
void get_load_channels(int K,int *load_channels,DiskReadMda &cluster) {
	for (int k=0; k<=K; k++) load_channels[k]=0;
	for (int i=0; i<cluster.N2(); i++) {
		int ch=(int)cluster.value(0,i);
		int k=(int)cluster.value(2,i);
		if ((1<=k)&&(k<=K)) load_channels[k]=ch;
	}
}
*/
