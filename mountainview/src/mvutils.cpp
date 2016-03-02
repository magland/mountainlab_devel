#include "mvutils.h"

#include <QCoreApplication>
#include "get_principal_components.h"
#include <math.h>

Mda compute_mean_waveform(DiskArrayModel *C) {
	Mda ret;
	if (!C->dim3()) return ret;
	int M=C->size(0);
	int T=C->size(1)/C->dim3();
	int NC=C->dim3();
	if (!NC) return ret;

	double sum[M*T];
	for (int ii=0; ii<M*T; ii++) sum[ii]=0;
	for (int c=0; c<NC; c++) {
		if ((c%100==0)||(c==NC-1)) {
			qApp->processEvents();
			//int pct=(int)(c*1.0/NC*100);
			//printf("Computing mean waveform...%d/%d (%d%%)\n",c,NC,pct);
		}
		int ii=0;
		Mda tmp=C->loadData(1,T*c,T*(c+1));
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				sum[ii]+=tmp.value(m,t);
				ii++;
			}
		}
	}
	ret.allocate(M,T);
	{
		int ii=0;
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				ret.setValue(sum[ii]/NC,m,t);
				ii++;
			}
		}
	}
	return ret;
}

Mda compute_mean_stdev_waveform(DiskArrayModel *C) {
	Mda ret;
	if (!C->dim3()) return ret;
	int M=C->size(0);
	int T=C->size(1)/C->dim3();
	int NC=C->dim3();
	if (!NC) return ret;

	double sum[M*T];
	double sumsqr[M*T];
	for (int ii=0; ii<M*T; ii++) {sum[ii]=0; sumsqr[ii]=0;}
	for (int c=0; c<NC; c++) {
		if ((c%100==0)||(c==NC-1)) {
			qApp->processEvents();
			//int pct=(int)(c*1.0/NC*100);
			//printf("Computing mean waveform...%d/%d (%d%%)\n",c,NC,pct);
		}
		int ii=0;
		Mda tmp=C->loadData(1,T*c,T*(c+1));
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				float val=tmp.value(m,t);
				sum[ii]+=val;
				sumsqr[ii]+=val*val;
				ii++;
			}
		}
	}
	ret.allocate(M,T*2);
	{
		int ii=0;
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				float mu=sum[ii]/NC;
				float sigma=sqrt(sumsqr[ii]/NC-mu*mu);
				float tmp=0;
				if (mu>0) {
					tmp=mu-sigma;
					if (tmp<0) tmp=0;
				}
				else {
					tmp=mu+sigma;
					if (tmp>0) tmp=0;
				}
				ret.setValue(mu,m,t);
				ret.setValue(tmp,m,T+t);
				ii++;
			}
		}
	}
	return ret;
}

Mda compute_features(DiskArrayModel *C) {
	Mda ret;
	if (!C->dim3()) return ret;
	int M=C->size(0);
	int T=C->size(1)/C->dim3();
	int NC=C->dim3();
	if (!NC) return ret;

	Mda X=C->loadData(1,0,T*NC);
	ret.allocate(3,NC);
	get_pca_features(M*T,NC,3,ret.dataPtr(),X.dataPtr());

	return ret;
}

Mda compute_features(const QList<DiskArrayModel *> &C) {
	Mda ret;
	if (C.isEmpty()) return ret;
	if (!C[0]->dim3()) return ret;
	int M=C[0]->size(0);
	int T=C[0]->size(1)/C[0]->dim3();
	int NC=0;
	for (int i=0; i<C.count(); i++) {
		NC+=C[i]->dim3();
	}
	if (!NC) return ret;

	Mda X; X.allocate(M,T,NC);
	int jj=0;
	for (int i=0; i<C.count(); i++) {
		Mda tmp=C[i]->loadData(1,0,T*C[i]->dim3());
		int ii=0;
		for (int aa=0; aa<T*C[i]->dim3(); aa++) {
			for (int m=0; m<M; m++) {
				X.setValue1(tmp.value1(ii),jj);
				ii++;
				jj++;
			}
		}
	}
	ret.allocate(3,NC);
	get_pca_features(M*T,NC,3,ret.dataPtr(),X.dataPtr());

	return ret;
}

