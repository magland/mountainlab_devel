#include "fit.h"
#include <stdio.h>
#include "mdaio.h"
#include <math.h>
#include <QTime>
#include <QList>
#include "mda.h"
#include "diskreadmda.h"

void sort_cluster_by_time(Mda &cluster);
void compute_scores(Mda &scores,Mda &X,Mda &templates,int num_events,int *times,int *labels,bool *score_update_needed,bool *coeffs_to_use);
Mda load_chunk(DiskReadMda &X,int tt1,int tt2);
int do_fit(bool *to_include,Mda &X,Mda &templates,int num_events,int *times,int *labels);

bool fit(const char *input_path,const char *templates_path,const char *cluster_in_path,const char *cluster_out_path) {
    printf("fit %s %s %s %s..\n",input_path,templates_path,cluster_in_path,cluster_out_path);

    DiskReadMda input;
    input.setPath(input_path);
    Mda cluster; cluster.read(cluster_in_path);
    Mda templates; templates.read(templates_path);

    sort_cluster_by_time(cluster);

    int num_events=cluster.N2();

    bool to_include[num_events];
    int times[num_events];
    int labels[num_events];
    for (int i=0; i<num_events; i++) {
        to_include[i]=0;
        times[i]=cluster.value(1,i)-1; //convert to zero-based indexing
        labels[i]=cluster.value(2,i);
    }

    int chunk_size=1000000;
    int num_chunks=(input.N2()+1)/chunk_size + 1;
    int overlap=10000;

    for (int chunk=0; chunk<num_chunks; chunk++) {
        int timepoint1=chunk*chunk_size;
        int timepoint2=timepoint1+chunk_size;
        if (timepoint2>input.N2()) timepoint2=input.N2();
        int tt1=timepoint1-overlap; if (tt1<0) tt1=0;
        int tt2=timepoint2+overlap; if (tt2>input.N2()) tt2=input.N2();
        Mda X=load_chunk(input,tt1,tt2);
        int num_events0=0;
        for (int j=0; j<num_events; j++) {
            if ((times[j]>=tt1)&&(times[j]<tt2)) num_events0++;
        }
        int times0[num_events0];
        int labels0[num_events0];
        bool to_include0[num_events];
        {
            int ii=0;
            for (int j=0; j<num_events; j++) {
                if ((times[j]>=tt1)&&(times[j]<tt2)) {
                    times0[ii]=times[j]-tt1;
                    labels0[ii]=labels[j];
                    ii++;
                }
            }
        }
        printf("chunk %d/%d %d events... ",chunk,num_chunks,num_events0);
        int num_passes=do_fit(to_include0,X,templates,num_events0,times0,labels0);
        int num_events_used=0;
        {
            int ii=0;
            for (int j=0; j<num_events; j++) {
                if ((times[j]>=tt1)&&(times[j]<tt2)) {
                    if (to_include0[ii]) num_events_used++;
                    if ((times[j]>=timepoint1)&&(times[j]<timepoint2)) {
                        to_include[j]=to_include0[ii];
                    }
                    ii++;
                }
            }
        }
        printf("%d events used in %d passes\n",num_events_used,num_passes);
    }

    int num_to_include=0;
    for (int j=0; j<num_events; j++) {
        if (to_include[j]) num_to_include++;
    }

    printf("Including %d / %d events.\n",num_to_include,num_events);
    Mda cluster_out; cluster_out.allocate(3,num_to_include);
    int i=0;
    for (int j=0; j<num_events; j++) {
        if (to_include[j]) {
            cluster_out.setValue(cluster.value(0,j),0,i);
            cluster_out.setValue(cluster.value(1,j),1,i);
            cluster_out.setValue(cluster.value(2,j),2,i);
            i++;
        }
    }

    cluster_out.write64(cluster_out_path);

    return true;
}

int do_fit(bool *to_include,Mda &X,Mda &templates,int num_events,int *times,int *labels) {
    int M=X.N1();
    int clip_size=templates.N2();
    int tt1=-(int)(clip_size/2);
    int tt2=tt1+clip_size-1;
    int K=templates.N3();

    for (int j=0; j<num_events; j++) {
        to_include[j]=false;
    }

	bool coeffs_to_use[M*clip_size*K];
	for (int k=0; k<K; k++) {
		float maxval=0;
		for (int t=0; t<clip_size; t++) {
			for (int m=0; m<M; m++) {
				float val=templates.value(m,t,k);
				if (fabs(val)>maxval) maxval=fabs(val);
			}
		}
		for (int t=0; t<clip_size; t++) {
			for (int m=0; m<M; m++) {
				float val=templates.value(m,t,k);
				coeffs_to_use[m+M*t+M*clip_size*k]=(fabs(val)>=maxval*0.2);
			}
		}
	}

    bool score_update_needed[num_events];
    for (int j=0; j<num_events; j++) score_update_needed[j]=true;

    int num_passes=0;
    Mda scores; scores.allocate(1,num_events);
    while (true) {
        num_passes++;
		compute_scores(scores,X,templates,num_events,times,labels,score_update_needed,coeffs_to_use);
        for (int j=0; j<num_events; j++) score_update_needed[j]=false;
        int num_changed=0;
		#pragma omp parallel for
        for (int j=0; j<num_events; j++) {
            if ((!to_include[j])&&(scores.value1(j)>0)) {
                bool best=true;
                int k=j-1;
                while ((k>=0)&&(fabs(times[k]-times[j])<=clip_size)) {
                    if (scores.value1(k)>scores.value1(j)) {
                        best=0;
                        break;
                    }
                    k--;
                }
                if (best) {
                    int k=j+1;
                    while ((k<num_events)&&(fabs(times[k]-times[j])<=clip_size)) {
                        if (scores.value1(k)>scores.value1(j)) {
                            best=0;
                            break;
                        }
                        k++;
                    }
                }
                if (best) {
                    if ((times[j]+tt1>=0)&&(times[j]+tt2<X.N2())) {
                        to_include[j]=true;
                        num_changed++;
                        for (int tt=tt1; tt<=tt2; tt++) {
                            for (int m=0; m<M; m++) {
                                float val1=X.value(m,times[j]+tt);
                                val1-=templates.value(m,tt-tt1,labels[j]-1);
                                X.setValue(val1,m,times[j]+tt);
                            }
                        }

                        {
                            int k=j;
                            while ((k>=0)&&(fabs(times[k]-times[j])<=clip_size)) {
                                score_update_needed[k]=true;
                                k--;
                            }
                        }
                        {
                            int k=j+1;
                            while ((k>=0)&&(fabs(times[k]-times[j])<=clip_size)) {
                                score_update_needed[k]=true;
                                k--;
                            }
                        }
                    }
                }
            }
        }
        if (num_changed==0) break;
    }
    return num_passes;
}

void compute_scores(Mda &scores,Mda &X,Mda &templates,int num_events,int *times,int *labels,bool *score_update_needed,bool *coeffs_to_use) {
    int M=X.N1();
    int clip_size=templates.N2();
    int tt1=-(int)(clip_size/2);
    int tt2=tt1+clip_size-1;

	#pragma omp parallel for
    for (int j=0; j<num_events; j++) {
        if (score_update_needed[j]) {
            if ((times[j]+tt1>=0)&&(times[j]+tt2<X.N2())) {
				int k=labels[j]-1;
				bool *to_use=&coeffs_to_use[M*clip_size*k];
                double sumsqr1=0;
                double sumsqr2=0;
                for (int tt=0; tt<clip_size; tt++) {
                    for (int m=0; m<M; m++) {
						if (to_use[m+M*tt]) {
							float X0=X.value(m,times[j]+tt1+tt);
							float T0=templates.value(m,tt,k);
							sumsqr1+=X0*X0;
							sumsqr2+=(X0-T0)*(X0-T0);
						}
                    }
                }
                scores.setValue1(sumsqr1-sumsqr2,j);
            }
            else scores.setValue1(0,j);
        }
    }
}

struct SortRec {
    double val;
    int index;
};
bool caseInsensitiveLessThan2(const SortRec &R1, const SortRec &R2)
{
  return R1.val<R2.val;
}


void sort_cluster_by_time(Mda &cluster) {
    QList<SortRec> list;
    for (int i=0; i<cluster.N2(); i++) {
        SortRec SR;
        SR.index=i;
        SR.val=cluster.value(1,i);
        list << SR;
    }
    qSort(list.begin(),list.end(),caseInsensitiveLessThan2);
    Mda cluster2; cluster2.allocate(cluster.N1(),cluster.N2());

    for (int i=0; i<list.count(); i++) {
        for (int jj=0; jj<cluster.N1(); jj++) {
            cluster2.setValue(cluster.value(jj,list[i].index),jj,i);
        }
    }

    cluster=cluster2;
}

Mda load_chunk(DiskReadMda &X,int tt1,int tt2) {
    Mda ret;
    int M=X.N1();
    int NN=tt2-tt1;
    ret.allocate(M,NN);
    for (int n=0; n<NN; n++) {
        for (int m=0; m<M; m++) {
            ret.setValue(X.value(m,n+tt1),m,n);
        }
    }
    return ret;
}
