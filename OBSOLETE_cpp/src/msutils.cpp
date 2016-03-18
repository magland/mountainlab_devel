#include "msutils.h"
#include <math.h>

double compute_min(const QList<double> &X) {
	double ret=X.value(0);
	for (int i=0; i<X.count(); i++) if (X[i]<ret) ret=X[i];
	return ret;
}

double compute_max(const QList<double> &X) {
	double ret=X.value(0);
	for (int i=0; i<X.count(); i++) if (X[i]>ret) ret=X[i];
	return ret;
}

int compute_max(const QList<int> &X) {
	int ret=X.value(0);
	for (int i=0; i<X.count(); i++) if (X[i]>ret) ret=X[i];
	return ret;
}

Mda extract_clips(DiskReadMda &X,const QList<double> &times,int clip_size) {
    QList<long> times2;
    for (int i=0; i<times.count(); i++) {
        times2 << times[i];
    }
    return extract_clips(X,times2,clip_size);
}

Mda extract_clips(DiskReadMda &X,const QList<long> &times,int clip_size) {
	int N=X.N2();
	int M=X.N1();
	int T=clip_size;
	int L=times.count();
	Mda clips; clips.allocate(M,T,L);
	double *clips_ptr=clips.dataPtr();
	int aaa=0;
	for (int i=0; i<L; i++) {
		int tt=times[i]-((T+1)/2 - 1);
		if ((tt>=0)&&(tt+T<=N)) {
			for (int t=0; t<T; t++) {
				for (int m=0; m<M; m++) {
					clips_ptr[aaa]=X.value1(m+M*(t+tt));
					//clips.setValue(,m,t,i);
					aaa++;
				}
			}
		}
		else aaa+=M*T; //important and tricky!!!
	}
	return clips;
}


Mda extract_clips(DiskReadMda &X,const QList<long> &times,const QList<int> &channels,int clip_size) {
	int N=X.N2();
	int M=channels.count();
	int T=clip_size;
	int L=times.count();
	Mda clips; clips.allocate(M,T,L);
	for (int i=0; i<L; i++) {
		int tt=times[i]-((T+1)/2 - 1);
		if ((tt>=0)&&(tt+T<=N)) {
			for (int t=0; t<T; t++) {
				for (int m=0; m<M; m++) {
					clips.setValue(X.value(channels[m],t+tt),m,t,i);
				}
			}
		}
	}
	return clips;
}

Mda compute_mean_clip(Mda &clips) {
	int M=clips.N1();
	int T=clips.N2();
	int L=clips.N3();
	Mda ret; ret.allocate(M,T);
	int aaa=0;
	for (int i=0; i<L; i++) {
		int bbb=0;
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				ret.setValue1(ret.value1(bbb)+clips.value1(aaa),bbb);
				aaa++;
				bbb++;
			}
		}
	}
	if (L) {
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				ret.setValue(ret.value(m,t)/L,m,t);
			}
		}
	}
	return ret;
}

double compute_mean(const QList<double> &X)
{
    double sum=0;
    for (int i=0; i<X.count(); i++) sum+=X[i];
    if (X.count()) sum/=X.count();
    return sum;
}

double compute_stdev(const QList<double> &X)
{
    double sumsqr=0;
    for (int i=0; i<X.count(); i++) sumsqr+=X[i]*X[i];
    double sum=0;
    for (int i=0; i<X.count(); i++) sum+=X[i];
    int ct=X.count();
    if (ct>=2) {
        return sqrt((sumsqr-sum*sum/ct)/(ct-1));
    }
    else return 0;
}
Mda grab_clips_subset(Mda &clips,const QList<int> &inds) {
    int M=clips.N1();
    int T=clips.N2();
    int LLL=inds.count();
    Mda ret; ret.allocate(M,T,LLL);
    for (int i=0; i<LLL; i++) {
        long aaa=i*M*T;
        long bbb=inds[i]*M*T;
        for (int k=0; k<M*T; k++) {
            ret.setValue1(clips.value1(bbb),aaa);
            aaa++; bbb++;
        }
    }
    return ret;
}

double compute_max(long N, double *X)
{
    if (N==0) return 0;
    double ret=X[0];
    for (long i=0; i<N; i++) {
        if (X[i]>ret) ret=X[i];
    }
    return ret;
}
