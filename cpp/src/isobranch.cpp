#include "isobranch.h"
#include "mda.h"
#include "isosplit.h"
#include "get_principal_components.h"
#include <QVector>

QVector<int> do_isobranch(Mda &clips, int min_cluster_size, int num_features, float isocut_threshold, int K_init);

bool isobranch(const char *clips_input_path, const char *labels_output_path, int min_cluster_size, int num_features, float isocut_threshold, int K_init)
{
    Mda clips;
    clips.read(clips_input_path);

    QVector<int> labels=do_isobranch(clips,min_cluster_size,num_features,isocut_threshold,K_init);
    Mda labels0;
    labels0.allocate(1,labels.count());
    for (int i=0; i<labels.count(); i++) {
        labels0.setValue(labels[i],0,i);
    }

    labels0.write(labels_output_path);

    return true;
}

int max_label(const QVector<int> &labels) {
    int ret=0;
    for (int i=0; i<labels.count(); i++) {
        if (labels[i]>ret) ret=labels[i];
    }
    return ret;
}

QList<int> get_label_inds(const QVector<int> &labels,int k) {
    QList<int> ret;
    for (int i=0; i<labels.count(); i++) {
        if (labels[i]==k) ret << i;
    }
    return ret;
}

Mda get_clips_subset(const Mda &clips,QList<int> &inds) {
    int M=clips.N1();
    int T=clips.N2();
    Mda ret;
    ret.allocate(M,T,inds.count());
    for (int i=0; i<inds.count(); i++) {
        for (int t=0; t<T; t++) {
            for (int m=0; m<M; m++) {
                ret.setValue(clips.value(m,t,inds[i]),m,t,i);
            }
        }
    }
    return ret;
}

QVector<float> compute_dists(Mda &clips) {
    int M=clips.N1();
    int T=clips.N2();
    int NC=clips.N3();

    int Tcenter=(int)((T+1)/2);

    QVector<float> ret;

    for (int i=0; i<NC; i++) {
        float maxval=0;
        for (int m=0; m<M; m++) {
            float val=fabs(clips.value(m,Tcenter,i));
            if (val>maxval) maxval=val;
        }
        ret << maxval;
    }

    return ret;

}

QVector<int> do_isobranch(Mda &clips, int min_cluster_size, int num_features, float isocut_threshold, int K_init) {
    int M=clips.N1();
    int T=clips.N2();
    int NC=clips.N3();

    QVector<int> all_ones;
    for (int i=0; i<NC; i++) all_ones << 1;

    if (NC<=min_cluster_size*2) {
        return all_ones;
    }

    Mda features;
    features.allocate(num_features,NC);
    get_pca_features(M*T,NC,num_features,features.dataPtr(),clips.dataPtr());

    QVector<int> labels=isosplit(features,isocut_threshold,K_init);
    int K=max_label(labels);
    printf("Found %d clusters. Size=%d\n",K,labels.count());

    if (K>1) {
        // There are multiple clusters, so let's divide and conquer
        QVector<int> labels_ret;
        for (int aa=0; aa<NC; aa++) labels_ret << 0;
        for (int k=0; k<K; k++) {
            QList<int> inds=get_label_inds(labels,k);
            int label_index=1;
            if (inds.count()>0) {
                Mda clips0=get_clips_subset(clips,inds);
                QVector<int> labels0=do_isobranch(clips0,min_cluster_size,num_features,isocut_threshold,K_init);
                int K0=max_label(labels0);
                for (int aa=0; aa<labels0.count(); aa++) {
                    if (labels0[aa]>0) {
                        labels_ret[inds[aa]]=labels0[aa]+label_index-1;
                    }
                }
                label_index+=K0;
            }
        }
        return labels_ret;
    }
    else {
        //we have only one cluster, so we need to remove some clips and repeat
        QVector<float> dists=compute_dists(clips);
        QVector<float> dists_sorted=dists; qSort(dists_sorted);
        float cutoff=dists_sorted[min_cluster_size];
        QList<int> inds;
        for (int i=0; i<NC; i++) {
            if (dists[i]<=cutoff) {
                inds << i;
            }
        }
        Mda clips0=get_clips_subset(clips,inds);
        QVector<int> labels0=do_isobranch(clips0,min_cluster_size,num_features,isocut_threshold,K_init);
        int K0=max_label(labels0);
        if (K0==1) {
            //there really is only one cluster! let's use it.
            return all_ones;
        }
        else {
            //there are multiple clusters in the subset. So let's set the things we excluded as zero.
            QVector<int> ret;
            for (int i=0; i<NC; i++) ret << 0;
            for (int i=0; i<inds.count(); i++) {
                ret[inds[i]]=labels[i];
            }
            return ret;
        }
    }
}
