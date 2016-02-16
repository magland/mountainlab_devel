#include "isobranch.h"
#include "mda.h"
#include "isosplit.h"
#include "get_principal_components.h"
#include <QVector>
#include <stdio.h>
#include <math.h>


QVector<int> do_isobranch(Mda &clips, const QList<float> &branch_thresholds, int num_features, float isocut_threshold, int K_init);

int max_label(const QVector<int> &labels) {
    int ret=0;
    for (int i=0; i<labels.count(); i++) {
        if (labels[i]>ret) ret=labels[i];
    }
    return ret;
}

int count_nonzero_labels(const QVector<int> &labels) {
    int ret=0;
    for (int i=0; i<labels.count(); i++) {
        if (labels[i]>0) ret++;
    }
    return ret;
}


bool isobranch(const char *clips_input_path, const char *labels_output_path, const QList<float> &branch_thresholds, int num_features, float isocut_threshold, int K_init)
{
    Mda clips;
    clips.read(clips_input_path);

    QVector<int> labels=do_isobranch(clips,branch_thresholds,num_features,isocut_threshold,K_init);
    Mda labels0;
    labels0.allocate(1,labels.count());
    for (int i=0; i<labels.count(); i++) {
        labels0.setValue(labels[i],0,i);
    }

    labels0.write(labels_output_path);

    printf("Using %d of %d events in %d clusters.\n",count_nonzero_labels(labels),labels.count(),max_label(labels));

    return true;
}


int compute_cluster_size(const QVector<int> &labels,int k) {
    int ret=0;
    for (int i=0; i<labels.count(); i++) {
        if (labels[i]==k) ret++;
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

QList<int> get_clips_above_threshold(const Mda &clips,float threshold) {
    int M=clips.N1();
    int T=clips.N2();
    int NC=clips.N3();

    int Tcenter=(int)((T+1)/2);

    QList<int> ret;

    for (int i=0; i<NC; i++) {
        float maxval=0;
        for (int m=0; m<M; m++) {
            float val=fabs(clips.value(m,Tcenter,i));
            if (val>maxval) maxval=val;
        }
        if (maxval>threshold) ret << i;
    }

    return ret;
}

struct Subclustering_node {
    int parent;
    QList<int> children;
};

struct Subclustering {
    QVector<int> labels;
    int K;
    QList<Subclustering_node> nodes;
};

QVector<int> compute_cluster_parents(int K,const QVector<int> &labels,int K_prev,const QVector<int> &labels_prev) {
    QVector<int> ret(K+1);
    if (K_prev==0) return ret;
    QVector<int> counts((K+1)*(K_prev+1));
    for (int ii=0; ii<labels.count(); ii++) {
        if ((labels[ii]>0)&&(labels_prev[ii]>0)) {
            counts[labels[ii]+(K+1)*labels_prev[ii]]++;
        }
    }
    for (int k=1; k<=K; k++) {
        int best_k_prev=1;
        int best_count=0;
        for (int k_prev=1; k_prev<=K_prev; k_prev++) {
            int ct=counts[k+(K+1)*k_prev];
            if (ct>best_count) {
                best_count=ct;
                best_k_prev=k_prev;
            }
        }
        ret[k]=best_k_prev;
    }
    return ret;
}

void create_tree_graph(QList<Subclustering> &subclusterings) {
    if (subclusterings.count()==0) return;
    for (int j=1; j<subclusterings.count(); j++) {
        Subclustering *SC=&subclusterings[j];
        Subclustering *SC_prev=&subclusterings[j-1];
        QVector<int> parents=compute_cluster_parents(SC->K,SC->labels,SC_prev->K,SC_prev->labels);
        for (int k=1; k<=SC->K; k++) {
            printf("C(%d,%d)->C(%d,%d)  %d<%d \n",j,k,j-1,parents[k],compute_cluster_size(SC->labels,k),compute_cluster_size(SC_prev->labels,parents[k]));
            SC->nodes[k].parent=parents[k];
            SC_prev->nodes[parents[k]].children << k;
        }
    }
}

bool is_nonbranching_node(const QList<Subclustering> &subclusterings,int j,int k) {
    if (subclusterings[j].nodes[k].children.count()==0) return true; //leaf node
    if (subclusterings[j].nodes[k].children.count()>1) return false; //more than one child means branching
    return is_nonbranching_node(subclusterings,j+1,subclusterings[j].nodes[k].children[0]); //check if only child is a nonbranching node
}

bool is_maximal_nonbranching_node(const QList<Subclustering> &subclusterings,int j,int k) {
    if (is_nonbranching_node(subclusterings,j,k)) {
        if (j==0) return true;
        int parent_k=subclusterings[j].nodes[k].parent;
        if (is_nonbranching_node(subclusterings,j-1,parent_k)) {
            return false; //not maximal
        }
        return true;
    }
    return false;
}

QVector<int> do_isobranch(Mda &clips,const QList<float> &branch_thresholds, int num_features, float isocut_threshold, int K_init) {
    int M=clips.N1();
    int T=clips.N2();
    int NC=clips.N3();

    if (branch_thresholds.isEmpty()) return QVector<int>(NC); //should not happen
    QList<float> branch_thresholds_next=branch_thresholds.mid(1);

    QList<int> inds0=get_clips_above_threshold(clips,branch_thresholds[0]);
    Mda clips0=get_clips_subset(clips,inds0);
    int NC0=clips0.N3();

    printf("Computing %d features of %d events at threshold %g...\n",num_features,NC0,branch_thresholds[0]);
    Mda features0; features0.allocate(num_features,NC0);
    get_pca_features(M*T,NC0,num_features,features0.dataPtr(),clips0.dataPtr());

    printf("ISO-SPLIT...\n");
    QVector<int> labels0=isosplit(features0,isocut_threshold,K_init);
    int K0=max_label(labels0);
    printf("Found %d clusters...\n",K0);

    QVector<int> ret(NC);
    int label_offset=0;
    for (int k=1; k<=K0; k++) {
        QList<int> inds1=get_label_inds(labels0,k);
        Mda clips1=get_clips_subset(clips0,inds1);
        QVector<int> labels1;
        if (!branch_thresholds_next.isEmpty())
            labels1=do_isobranch(clips1,branch_thresholds_next,num_features,isocut_threshold,K_init);
        else
            labels1=QVector<int>(inds1.count(),1); //all ones
        int K1=max_label(labels1);
        for (int i=0; i<inds1.count(); i++) {
            if (labels1[i]>0) {
                ret[inds0[inds1[i]]]=label_offset+labels1[i];
            }
        }
        label_offset+=K1;
    }
    if (max_label(ret)<=1) {
        //only one branch exists! That means we should just use this cluster as a whole
        ret=QVector<int>(NC,1);
    }

    return ret;
}

QVector<int> do_isobranch_old(Mda &clips, const QList<float> &branch_thresholds, int num_features, float isocut_threshold, int K_init) {
    int M=clips.N1();
    int T=clips.N2();
    int NC=clips.N3();

    QList<Subclustering> subclusterings;

    for (int j=0; j<branch_thresholds.count(); j++) {
        printf("Branch threshold = %f...\n",branch_thresholds[j]);
        Subclustering SC0;
        QList<int> inds=get_clips_above_threshold(clips,branch_thresholds[j]);
        Mda clips0=get_clips_subset(clips,inds);
        int NC0=clips0.N3();
        Mda features0; features0.allocate(num_features,NC0);
        printf("Computing %d features of %d events...\n",num_features,NC0);
        get_pca_features(M*T,NC0,num_features,features0.dataPtr(),clips0.dataPtr());
        printf("ISO-SPLIT...\n");
        QVector<int> labels0=isosplit(features0,isocut_threshold,K_init);
        SC0.labels=QVector<int>(NC);
        for (int ii=0; ii<inds.count(); ii++) {
            SC0.labels[inds[ii]]=labels0[ii];
        }
        SC0.K=max_label(SC0.labels);
        printf("Found %d clusters...\n",SC0.K);
        Subclustering_node dummy_node; dummy_node.parent=0;
        for (int k=0; k<=SC0.K; k++) SC0.nodes << dummy_node;
        subclusterings << SC0;
    }

    printf("Creating tree graph...\n");
    create_tree_graph(subclusterings);

    printf("Finalizing...\n");
    QVector<int> labels(NC);
    int label_number=1;
    for (int j=0; j<subclusterings.count(); j++) {
        for (int k=1; k<subclusterings[j].K; k++) {
            if (is_maximal_nonbranching_node(subclusterings,j,k)) {
                QList<int> inds=get_label_inds(subclusterings[j].labels,k);
                printf("Using cluster C(%d,%d) of size %d\n",j,k,inds.count());
                for (int i=0; i<inds.count(); i++) {
                    labels[inds[i]]=label_number;
                }
                label_number++;
            }
        }
    }
    printf("Using %d clusters and %d/%d events.\n",max_label(labels),count_nonzero_labels(labels),labels.count());

    return labels;
}


