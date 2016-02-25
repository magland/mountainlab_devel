#include "isosplit2.h"
#include <QSet>
#include <math.h>

QList<int> do_kmeans(Mda &X,int K);

struct AttemptedComparisons {
    QList<double> centers1,centers2;
    QList<int> counts1,counts2;
};

QList<int> find_inds(const QList<int> &labels,int k) {
    QList<int> ret;
    for (int i=0; i<labels.count(); i++) {
        if (labels[i]==k) ret << i;
    }
    return ret;
}

void geometric_median(int M,int N,double *ret,double *X) {
    int num_iterations=10;
    if (N==0) return;
    if (N==1) {
        for (int m=0; m<M; m++) ret[m]=X[m];
        return;
    }
    double *weights=(double *)malloc(sizeof(double)*N);
    for (int i=0; i<N; i++) weights[i]=1;
    for (int it=1; it<=num_iterations; it++) {
        double sum_weights=0;
        for (int i=0; i<N; i++) sum_weights+=weights[i];
        if (sum_weights) {
            for (int i=0; i<N; i++) weights[i]/=sum_weights;
        }
        for (int m=0; m<M; m++) ret[m]=0;
        int aa=0;
        for (int n=0; n<N; n++) {
            for (int m=0; m<M; m++) {
                ret[m]+=X[aa]*weights[n];
                aa++;
            }
        }
        aa=0;
        for (int n=0; n<N; n++) {
            double sumsqr_diff=0;
            for (int m=0; m<M; m++) {
                double val=X[aa]-ret[m];
                sumsqr_diff=val*val;
                aa++;
            }
            if (sumsqr_diff!=0) {
                weights[n]=1/sqrt(sumsqr_diff);
            }
            else weights[n]=0;
        }
    }
}

QList<double> compute_center(Mda &X,const QList<int> &inds) {
    int M=X.N1();
    int NN=inds.count();
    if (NN==0) {
        QList<double> ret; for (int i=0; i<M; i++) ret << 0;
        return ret;
    }
    double *XX=(double *)malloc(sizeof(double)*NN);
    int aa=0;
    for (int n=0; n<NN; n++) {
        for (int m=0; m<M; m++) {
            XX[aa]=X.value(m,inds[n]);
            aa++;
        }
    }
    double *result=(double *)malloc(sizeof(double)*M);
    geometric_median(M,NN,result,XX);
    QList<double> ret; for (int m=0; m<M; m++) ret << result[m];
    free(result);
    free(XX);
    return ret;
}

Mda compute_centers(Mda &X,const QList<int> &labels,int K) {
    int M=X.N1();
    //int N=X.N2();
    Mda ret;
    ret.allocate(M,K);
    for (int k=0; k<K; k++) {
        QList<int> inds=find_inds(labels,k);
        QList<double> ctr=compute_center(X,inds);
        for (int m=0; m<M; m++) ret.setValue(ctr[m],m,k);
    }
    return ret;
}

double distance_between_vectors(int M,double *v1,double *v2) {
    double sumsqr=0;
    for (int i=0; i<M; i++) {
        double val=v1[i]-v2[i];
        sumsqr+=val*val;
    }
    return sqrt(sumsqr);
}


bool was_already_attempted(int M,AttemptedComparisons &attempted_comparisons,double *center1,double *center2,int count1,int count2,double repeat_tolerance) {
    double tol=repeat_tolerance;
    for (int i=0; i<attempted_comparisons.counts1.count(); i++) {
        double diff_count1=fabs(attempted_comparisons.counts1[i]-count1);
        double avg_count1=(attempted_comparisons.counts1[i]+count1)/2;
        if (diff_count1<=tol*avg_count1) {
            double diff_count2=fabs(attempted_comparisons.counts2[i]-count2);
            double avg_count2=(attempted_comparisons.counts2[i]+count2)/2;
            if (diff_count2<=tol*avg_count2) {
                double C1[M]; for (int m=0; m<M; m++) C1[m]=attempted_comparisons.centers1[i*M+m];
                double C2[M]; for (int m=0; m<M; m++) C2[m]=attempted_comparisons.centers2[i*M+m];
                double dist0=distance_between_vectors(M,C1,C2);
                double dist1=distance_between_vectors(M,C1,center1);
                double dist2=distance_between_vectors(M,C2,center2);
                if (dist0>0) {
                    double frac1=dist1/dist0;
                    double frac2=dist2/dist0;
                    if ((frac1<=tol*1/sqrt(count1))&&(frac2<=tol*1/sqrt(count2))) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

void find_next_comparison(int M,int K,int &k1,int &k2,bool *active_labels,double *Cptr,int *counts,AttemptedComparisons &attempted_comparisons,double repeat_tolerance) {
    QList<int> active_inds;
    for (int k=0; k<K; k++) if (active_labels[k]) active_inds << k;
    if (active_inds.count()==0) {
        k1=-1; k2=-1;
        return;
    }
    int Nactive=active_inds.count();
    double dists[Nactive][Nactive];
    for (int a=0; a<Nactive; a++) {
        for (int b=0; b<Nactive; b++) {
            dists[a][b]=distance_between_vectors(M,&Cptr[active_inds[a]*M],&Cptr[active_inds[b]*M]);
        }
        dists[a][a]=-1; //don't use it
    }
    while (true) {
        int best_a=-1,best_b=-1;
        double best_dist=-1;
        for (int a=0; a<Nactive; a++) {
            for (int b=0; b<Nactive; b++) {
                double val=dists[a][b];
                if (val>=0) {
                    if ((best_dist<0)||(val<best_dist)) {
                        best_a=a;
                        best_b=b;
                        best_dist=val;
                    }
                }
            }
        }
        k1=active_inds[best_a];
        k2=active_inds[best_b];
        if ((counts[k1]>0)&&(counts[k2]>0)) { //just to make sure (zero was actually happening some times, but not sure why)
            if (!was_already_attempted(M,attempted_comparisons,&Cptr[k1*M],&Cptr[k2*M],counts[k1],counts[k2],repeat_tolerance)) {
                //hurray!
                return;
            }
        }
        dists[k1][k2]=-1;
        dists[k2][k1]=-1;
    }
    k1=-1; k2=-1;
}

QList<int> test_redistribute(bool &do_merge,Mda &X1,Mda &X2,double isocut_threshold) {
//    if opts.whiten_at_each_comparison
//        [X1,X2,V]=whiten_two_clusters(X1,X2);
//    else
//        V=compute_cluster_center(X2)-compute_cluster_center(X1);
//    end;

//    if (sum(V.^2)==0)
//        warning('isosplit: vector V is null.');
//    else
//        V=V/sqrt(sum(V.^2));
//    end;
//    XX=V'*cat(2,X1,X2); %Project onto the line connecting the centroids
//    N=length(XX);
//    if (N<=5) %avoid a crash - 2/22/2016 jfm
//        do_merge=1;
//        labels=ones(1,N);
//        return;
//    end;
//    XXs=sort(XX);
//    cutpoint=isocut(XXs,opts.isocut_threshold); %This is the core procedure -- split based on isotonic regression
//    if (cutpoint~=0)
//        %It was a statistically significant split -- so let's redistribute!
//        ii1=find(XX<=cutpoint);
//        ii2=find(XX>cutpoint);
//        do_merge=0;
//    else
//        ii1=1:N;
//        ii2=[];
//        do_merge=1;
//    end;
//    labels=zeros(1,N);
//    labels(ii1)=1;
//    labels(ii2)=2;
}

QList<int> test_redistribute(bool &do_merge,Mda &X,const QList<int> &inds1,const QList<int> &inds2,double isocut_threshold) {
    int M=X.N1();
    Mda X1;
    X1.allocate(M,inds1.count());
    for (int i=0; i<inds1.count(); i++) {
        for (int m=0; m<M; m++) {
            X1.setValue(X.value(m,inds1[i]),m,i);
        }
    }
    Mda X2;
    X2.allocate(M,inds2.count());
    for (int i=0; i<inds2.count(); i++) {
        for (int m=0; m<M; m++) {
            X2.setValue(X.value(m,inds2[i]),m,i);
        }
    }
    return test_redistribute(do_merge,X1,X2,isocut_threshold);
}

int compute_max_00(const QList<int> &X) {
    int ret=X.value(0);
    for (int i=0; i<X.count(); i++) if (X[i]>ret) ret=X[i];
    return ret;
}


QList<int> isosplit2(Mda &X, float isocut_threshold, int K_init)
{
    double repeat_tolerance=0.2;

    int M=X.N1();
    int N=X.N2();
    QList<int> labels=do_kmeans(X,K_init);

    bool active_labels[K_init];
    for (int ii=0; ii<K_init; ii++) active_labels[ii]=true;
    Mda centers=compute_centers(X,labels,K_init); //M x K_init
    int counts[K_init]; for (int ii=0; ii<K_init; ii++) counts[ii]=0;
    for (int i=0; i<N; i++) counts[labels[i]]++;
    double *Cptr=centers.dataPtr();

    AttemptedComparisons attempted_comparisons;

    int num_iterations=0;
    int max_iterations=1000;
    while ((true)&&(num_iterations<max_iterations)) {
        num_iterations++;
        QList<int> old_labels=labels;
        int k1,k2;
        find_next_comparison(M,K_init,k1,k2,active_labels,Cptr,counts,attempted_comparisons,repeat_tolerance);
        if (k1<0) break;

        QList<int> inds1=find_inds(labels,k1);
        QList<int> inds2=find_inds(labels,k2);
        QList<int> inds12=inds1; inds12.append(inds2);
        for (int m=0; m<M; m++) {
            attempted_comparisons.centers1 << Cptr[m+k1*M];
            attempted_comparisons.centers2 << Cptr[m+k2*M];
        }
        attempted_comparisons.counts1 << inds1.count();
        attempted_comparisons.counts2 << inds2.count();
        for (int m=0; m<M; m++) {
            attempted_comparisons.centers2 << Cptr[m+k1*M];
            attempted_comparisons.centers1 << Cptr[m+k2*M];
        }
        attempted_comparisons.counts2 << inds1.count();
        attempted_comparisons.counts1 << inds2.count();

        bool do_merge;
        QList<int> labels0=test_redistribute(do_merge,X,inds1,inds2,isocut_threshold);
        int max_label=compute_max_00(labels0);
        if ((do_merge)||(max_label==1)) {
            for (int i=0; i<N; i++) {
                if (labels[i]==k2) labels[i]=k1;
            }
            QList<double> ctr=compute_center(X,inds12);
            for (int m=0; m<M; m++) {
                centers.setValue(ctr[m],m,k1);
            }
            counts[k1]=inds12.count();
            counts[k2]=0;
            active_labels[k2]=false;

        }
        else {
            QList<int> indsA0=find_inds(labels0,1);
            QList<int> indsB0=find_inds(labels0,2);
            QList<int> indsA,indsB;
            for (int i=0; i<indsA0.count(); i++) indsA << inds12[indsA0[i]];
            for (int i=0; i<indsB0.count(); i++) indsB << inds12[indsB0[i]];
            for (int i=0; i<indsA.count(); i++) {
                labels[indsA[i]]=k1;
            }
            for (int i=0; i<indsB.count(); i++) {
                labels[indsB[i]]=k2;
            }
            {
                QList<double> ctr=compute_center(X,indsA);
                for (int m=0; m<M; m++) {
                    centers.setValue(ctr[m],m,k1);
                }
            }
            {
                QList<double> ctr=compute_center(X,indsB);
                for (int m=0; m<M; m++) {
                    centers.setValue(ctr[m],m,k2);
                }
            }
            counts[k1]=indsA.count();
            counts[k2]=indsB.count();
        }

    }

    int labels_map[K_init]; for (int k=0; k<K_init; k++) labels_map[k]=0;
    int kk=1;
    for (int j=0; j<K_init; j++) {
        if ((active_labels[j])&&(counts[j]>0)) {
            labels_map[j]=kk; kk++;
        }
    }
    QList<int> ret;
    for (int i=0; i<N; i++) {
        ret << labels_map[labels[i]];
    }
    return ret;
}

//choose K distinct (sorted) integers between 0 and N-1. If K>N then it will repeat the last integer a suitable number of times
QList<int> choose_random_indices(int N,int K) {;
    QList<int> ret;
    if (K>=N) {
        for (int i=0; i<N; i++) ret << i;
        while (ret.count()<K) ret << N-1;
        return ret;
    }
    QSet<int> theset;
    while (theset.count()<K) {
        int ind=(qrand()%N);
        theset.insert(ind);
    }
    ret=theset.toList();
    qSort(ret);
    return ret;
}

//do k-means with K clusters -- X is MxN representing N points in M-dimensional space. Returns a labels vector of size N.
QList<int> do_kmeans(Mda &X,int K) {
    int M=X.N1();
    int N=X.N2();
    double *Xptr=X.dataPtr();
    Mda centroids_mda; centroids_mda.allocate(M,K); double *centroids=centroids_mda.dataPtr();
    QList<int> labels; for (int i=0; i<N; i++) labels << -1;
    int *counts=(int *)malloc(sizeof(int)*K);

    //initialize the centroids
    QList<int> initial=choose_random_indices(N,K);
    for (int j=0; j<K; j++) {
        int ind=initial[j];
        int jj=ind*M;
        int ii=j*M;
        for (int m=0; m<M; m++) {
            centroids[m+ii]=Xptr[m+jj];
        }
    }

    bool something_changed=true;
    while (something_changed) {
        something_changed=false;
        //Assign the labels
        for (int n=0; n<N; n++) {
            int jj=n*M;
            double best_distsqr=0;
            int best_k=0;
            for (int k=0; k<K; k++) {
                int ii=k*M;
                double tmp=0;
                for (int m=0; m<M; m++) {
                    tmp+=(centroids[m+ii]-Xptr[m+jj])*(centroids[m+ii]-Xptr[m+jj]);
                }
                if ((k==0)||(tmp<best_distsqr)) {
                    best_distsqr=tmp;
                    best_k=k;
                }
            }
            if (labels[n]!=best_k) {
                labels[n]=best_k;
                something_changed=true;
            }
        }

        if (something_changed) {
            //Compute the centroids
            for (int k=0; k<K; k++) {
                int ii=k*M;
                for (int m=0; m<M; m++) {
                    centroids[m+ii]=0;
                }
                counts[k]=0;
            }
            for (int n=0; n<N; n++) {
                int jj=n*M;
                int k=labels[n];
                int ii=k*M;
                for (int m=0; m<M; m++) {
                    centroids[m+ii]+=Xptr[m+jj];
                }
                counts[k]++;
            }
            for (int k=0; k<K; k++) {
                int ii=k*M;
                if (counts[k]) {
                    for (int m=0; m<M; m++)
                        centroids[m+ii]/=counts[k];
                }
            }
        }
    }

    free(counts);

    return labels;
}
