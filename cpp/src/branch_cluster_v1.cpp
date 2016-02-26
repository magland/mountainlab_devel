#include "branch_cluster_v1.h"
#include "diskreadmda.h"
#include <stdio.h>
#include "get_principal_components.h"
#include <math.h>
#include "isosplit2.h"

Mda extract_clips_00(DiskReadMda &X,const QList<long> &times,const QList<int> &channels,int clip_size);
QList<int> do_branch_cluster(Mda &clips,const Branch_Cluster_Opts &opts,QList<int> &base_inds);
QList<double> compute_peaks(Mda &clips,int ch);
double compute_min(const QList<double> &X);
double compute_max(const QList<double> &X);
int compute_max(const QList<int> &X);
QList<int> consolidate_labels(DiskReadMda &X,const QList<long> &times,const QList<int> &labels,int ch,int clip_size);

bool branch_cluster_v1(const char *raw_path, const char *detect_path, const char *adjacency_matrix_path, const char *output_firings_path, const Branch_Cluster_Opts &opts)
{
    DiskReadMda X; X.setPath(raw_path);
    int M=X.N1();

    DiskReadMda detect; detect.setPath(detect_path);
    int L=detect.N2();

    Mda AM;
    if ((adjacency_matrix_path)&&(strlen(adjacency_matrix_path)>0)) {
        AM.read(adjacency_matrix_path);
    }
    else {
        AM.allocate(M,M);
        for (int i=0; i<M; i++) {
            for (int j=0; j<M; j++) {
                AM.setValue(1,i,j);
            }
        }
    }

    if ((AM.N1()!=M)||(AM.N2()!=M)) {
        printf("Error: incompatible dimensions between AM and X.\n");
        return false;
    }

    Mda firings0; firings0.allocate(5,L); //L is the max it could be

    int jjjj=0;
    int k_offset=0;
    for (int m=0; m<M; m++) {
        printf("Channel %d/%d...\n",m+1,M);
        QList<int> neighborhood; neighborhood << m;
        for (int a=0; a<M; a++) if ((AM.value(m,a))&&(a!=m)) neighborhood << a;
        QList<long> times;
        for (int i=0; i<L; i++) {
            if (detect.value(0,i)==(m+1)) {
                times << (long)detect.value(1,i) - 1; //convert to 0-based indexing
            }
        }
        Mda clips=extract_clips_00(X,times,neighborhood,opts.clip_size);
        QList<int> base_inds;
        QList<int> labels=do_branch_cluster(clips,opts,base_inds);
        printf("\n");
        labels=consolidate_labels(X,times,labels,m,opts.clip_size);
        QList<double> peaks=compute_peaks(clips,0);

        for (int i=0; i<times.count(); i++) {
            if (labels[i]) {
                firings0.setValue(m+1,0,jjjj); //channel
                firings0.setValue(times[i]+1,1,jjjj); //times //convert back to 1-based indexing
                firings0.setValue(labels[i]+k_offset,2,jjjj); //labels
                firings0.setValue(peaks[i],3,jjjj); //peaks
                jjjj++;
            }
        }
        k_offset+=compute_max(labels);
        printf("\n");
    }

    int L_true=jjjj;
    Mda firings; firings.allocate(firings0.N1(),L_true);
    for (int i=0; i<L_true; i++) {
        for (int j=0; j<firings0.N1(); j++) {
            firings.setValue(firings0.value(j,i),j,i);
        }
    }

    firings.write(output_firings_path);

    return true;
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

QList<int> consolidate_labels(DiskReadMda &X,const QList<long> &times,const QList<int> &labels,int ch,int clip_size) {
    int M=X.N1();
    int T=clip_size;
    int K=compute_max(labels);
    QList<int> all_channels;
    for (int m=0; m<M; m++) all_channels << m;
    int label_mapping[K+1];
    label_mapping[0]=0;
    int kk=1;
    for (int k=1; k<=K; k++) {
        QList<long> times_k;
        for (int i=0; i<times.count(); i++) {
            if (labels[i]==k) times_k << times[i];
        }
        Mda clips_k=extract_clips_00(X,times_k,all_channels,clip_size);
        Mda template_k=compute_mean_clip(clips_k);
        QList<double> energies;
        for (int m=0; m<M; m++) energies << 0;
        for (int t=0; t<T; t++) {
            for (int m=0; m<M; m++) {
                double val=template_k.value(m,t);
                energies[m]+=val*val;
            }
        }
        double max_energy=compute_max(energies);
        if (energies[ch]>=max_energy*0.9) {
            label_mapping[k]=kk;
            kk++;
        }
        else label_mapping[k]=0;
    }
    QList<int> ret;
    for (int i=0; i<labels.count(); i++) {
        ret << label_mapping[labels[i]];
    }
    printf("Using %d of %d clusters.\n",compute_max(ret),K);
    return ret;
}

QList<double> compute_peaks(Mda &clips,int ch) {
    int T=clips.N2();
    int L=clips.N3();
    int t0=(T+1)/2 - 1;
    QList<double> ret;
    for (int i=0; i<L; i++) {
        ret << clips.value(ch,t0,i);
    }
    return ret;
}

QList<double> compute_abs_peaks(Mda &clips,int ch) {
    int T=clips.N2();
    int L=clips.N3();
    int t0=(T+1)/2 - 1;
    QList<double> ret;
    for (int i=0; i<L; i++) {
        ret << fabs(clips.value(ch,t0,i));
    }
    return ret;
}

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

QList<int> find_peaks_below_threshold(QList<double> &peaks,double threshold) {
    QList<int> ret;
    for (int i=0; i<peaks.count(); i++) {
        if (peaks[i]<threshold) ret << i;
    }
    return ret;
}

QList<int> find_peaks_above_threshold(QList<double> &peaks,double threshold) {
    QList<int> ret;
    for (int i=0; i<peaks.count(); i++) {
        if (peaks[i]>=threshold) ret << i;
    }
    return ret;
}

void normalize_features(Mda &F) {
    int M=F.N1();
    int N=F.N2();
    QList<double> norms;
    int aa=0;
    for (int i=0; i<N; i++) {
        double sumsqr=0;
        for (int j=0; j<M; j++) {
            double val=F.value1(aa);
            sumsqr+=val*val;
            aa++;
        }
        norms << sqrt(sumsqr);
    }
    aa=0;
    for (int i=0; i<N; i++) {
        double factor=1;
        if (norms[i]) factor=1/norms[i];
        for (int j=0; j<M; j++) {
            double val=F.value1(aa);
            F.setValue1(val*factor,aa);
            aa++;
        }
    }
}

QList<int> do_cluster_with_normalized_features(Mda &clips,const Branch_Cluster_Opts &opts) {
    int M=clips.N1();
    int T=clips.N2();
    int L=clips.N3();
    int nF=opts.num_features;
    Mda FF; FF.allocate(nF,L);
    get_pca_features(M*T,L,nF,FF.dataPtr(),clips.dataPtr());
    normalize_features(FF);
    return isosplit2(FF);
}

QList<int> do_cluster_without_normalized_features(Mda &clips,const Branch_Cluster_Opts &opts) {
    int M=clips.N1();
    int T=clips.N2();
    int L=clips.N3();
    int nF=opts.num_features;
    Mda FF; FF.allocate(nF,L);
    get_pca_features(M*T,L,nF,FF.dataPtr(),clips.dataPtr());
    //normalize_features(FF);
    return isosplit2(FF);
}

QList<int> do_branch_cluster(Mda &clips,const Branch_Cluster_Opts &opts,QList<int> &base_inds) {
    int L=clips.N3();
    QList<double> peaks=compute_peaks(clips,0);
    QList<double> abs_peaks=compute_abs_peaks(clips,0);
    base_inds.clear();

    //In the case we have both positive and negative peaks, just split into two tasks!
    double min0=compute_min(peaks);
    double max0=compute_max(peaks);
    if ((min0<0)&&(max0>=0)) {
        QList<int> inds_neg,inds_pos;
        for (int i=0; i<L; i++) {
            if (peaks[i]<0) inds_neg << i;
            else inds_pos << i;
        }
        Mda clips_neg=grab_clips_subset(clips,inds_neg);
        Mda clips_pos=grab_clips_subset(clips,inds_pos);
        QList<int> base_inds_neg,base_inds_pos;
        printf("==NEGATIVES (%d)==\n",inds_neg.count());
        QList<int> labels_neg=do_branch_cluster(clips_neg,opts,base_inds_neg);
        printf("\n==POSITIVES (%d)==\n",inds_pos.count());
        QList<int> labels_pos=do_branch_cluster(clips_pos,opts,base_inds_pos);
        int K_neg=compute_max(labels_neg);
        QList<int> labels; for (int i=0; i<L; i++) labels << 0;
        for (int i=0; i<inds_neg.count(); i++) {
            labels[inds_neg[i]]=labels_neg[i];
        }
        for (int i=0; i<inds_pos.count(); i++) {
            if (labels_pos[i]) labels[inds_pos[i]]=labels_pos[i]+K_neg;
            else labels[inds_pos[i]]=0;
        }
        return labels;
    }

    //QList<int> labels0=do_cluster_with_normalized_features(clips,opts);
    QList<int> labels0=do_cluster_without_normalized_features(clips,opts);
    int K0=compute_max(labels0);
    if (K0>1) {
        printf("(K=%d:",K0);
        QList<int> labels; for (int i=0; i<L; i++) labels << 0;
        int kk_offset=0;
        for (int k=1; k<=K0; k++) {
            QList<int> inds_k;
            for (int a=0; a<L; a++) {
                if (labels0[a]==k) inds_k << a;
            }
            printf(" %d",inds_k.count());
            Mda clips_k=grab_clips_subset(clips,inds_k);
            QList<int> base_inds_k;
            QList<int> labels_k=do_branch_cluster(clips_k,opts,base_inds_k);
            for (int a=0; a<inds_k.count(); a++) {
                labels[inds_k[a]]=labels_k[a]+kk_offset;
            }
            kk_offset+=compute_max(labels_k);
        }
        printf(")");
        return labels;
    }
    else {
        double abs_peak_threshold=0;
        double max_abs_peak=compute_max(abs_peaks);
        while (true) {
            QList<int> inds_below=find_peaks_below_threshold(abs_peaks,abs_peak_threshold);
            if ((inds_below.count()>=opts.min_section_count)&&(L-inds_below.count()>=opts.min_section_count)) {
                break;
            }
            if (abs_peak_threshold>max_abs_peak) {
                break;
            }
            abs_peak_threshold+=opts.section_increment;
        }
        if (abs_peak_threshold>max_abs_peak) {
            QList<int> labels; for (int i=0; i<L; i++) labels << 1;
            for (int i=0; i<L; i++) base_inds << i;
            return labels;
        }
        else {
            QList<int> labels; for (int i=0; i<L; i++) labels << 1;
            QList<int> inds_below=find_peaks_below_threshold(abs_peaks,abs_peak_threshold);
            QList<int> inds_above=find_peaks_above_threshold(abs_peaks,abs_peak_threshold);
            Mda clips_above=grab_clips_subset(clips,inds_above);
            QList<int> base_inds_0;
            QList<int> labels_0=do_branch_cluster(clips_above,opts,base_inds_0);
            int k_offset=0;
            if (base_inds_0.count()==0) k_offset=1;
            for (int i=0; i<inds_above.count(); i++) {
                labels[inds_above[i]]=labels_0[i]+k_offset;
            }
            for (int i=0; i<inds_below.count(); i++) {
                base_inds << inds_below[i];
            }
            for (int i=0; i<base_inds_0.count(); i++) {
                base_inds << inds_above[base_inds_0[i]];
            }
            return labels;
        }
    }
}

Mda extract_clips_00(DiskReadMda &X,const QList<long> &times,const QList<int> &channels,int clip_size) {
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
