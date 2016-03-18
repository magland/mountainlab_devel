#include "remove_noise_subclusters.h"
#include "diskreadmda.h"
#include <math.h>
#include "get_pca_features.h"
#include <QTime>
#include <stdio.h>
#include <QDebug>
#include "msmisc.h"
#include "extract_clips.h"

struct Subcluster {
	QList<int> inds;
	double detectability_score;
	double peak;
	int label;
};

QList<int> find_label_inds(const QList<int> &labels,int k) {
	QList<int> ret;
	for (int i=0; i<labels.count(); i++) {
		if (labels[i]==k) ret << i;
	}
	return ret;
}

Mda get_subclips(Mda &clips,const QList<int> &inds) {
	int M=clips.N1();
	int T=clips.N2();
	int L2=inds.count();
	Mda ret; ret.allocate(M,T,L2);
	for (int i=0; i<L2; i++) {
		int aaa=M*T*i;
		int bbb=M*T*inds[i];
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				ret.set(clips.get(bbb),aaa);
				aaa++;
				bbb++;
			}
		}
	}
	return ret;
}

QList<Shell> define_shells(const QList<double> &peaks,const Define_Shells_Opts &opts) {
	QList<Shell> shells;
	if (peaks.count()==0) return shells;

	//negatives and positives
	if (compute_min(peaks)<0) {
		QList<int> inds_neg,inds_pos;
		for (int i=0; i<peaks.count(); i++) {
			if (peaks[i]<0) inds_neg << i;
			else inds_pos << i;
		}
		QList<double> neg_peaks_neg; for (int i=0; i<inds_neg.count(); i++) neg_peaks_neg << -peaks[inds_neg[i]];
		QList<double> peaks_pos; for (int i=0; i<inds_pos.count(); i++) peaks_pos << peaks[inds_pos[i]];
		QList<Shell> shells_neg=define_shells(neg_peaks_neg,opts);
		QList<Shell> shells_pos=define_shells(peaks_pos,opts);
		for (int i=shells_neg.count()-1; i>=0; i--) {
			Shell SS=shells_neg[i];
			QList<int> new_inds; for (int i=0; i<SS.inds.count(); i++) new_inds << inds_neg[SS.inds[i]];
			SS.inds=new_inds;
			shells << SS;
		}
		for (int i=0; i<shells_pos.count(); i++) {
			Shell SS=shells_pos[i];
			QList<int> new_inds; for (int i=0; i<SS.inds.count(); i++) new_inds << inds_pos[SS.inds[i]];
			SS.inds=new_inds;
			shells << SS;
		}
		return shells;
	}

	//only positives
	double max_peak=compute_max(peaks);
	double peak_lower=0;
	double peak_upper=peak_lower+opts.shell_increment;
	while (peak_lower<=max_peak) {
		QList<int> inds1;
		QList<int> inds2;
		for (int i=0; i<peaks.count(); i++) {
			if ((peak_lower<=peaks[i])&&(peaks[i]<peak_upper)) {
				inds1 << i;
			}
			if (peaks[i]>=peak_upper) {
				inds2 << i;
			}
		}
		int ct1=inds1.count();
		int ct2=inds2.count();
		if (peak_upper>max_peak) {
			Shell SSS; SSS.inds=inds1;
			peak_lower=peak_upper;
			shells << SSS;
		}
		else if ((ct1>=opts.min_shell_size)&&(ct2>=opts.min_shell_size)) {
			Shell SSS; SSS.inds=inds1;
			peak_lower=peak_upper;
			peak_upper=peak_lower+opts.shell_increment;
			shells << SSS;
		}
		else {
			peak_upper=peak_upper+opts.shell_increment;
		}
	}
	return shells;

}

QList<double> randsample_with_replacement(long N,long K) {
	QList<double> ret;
	for (int i=0; i<K; i++) {
		ret << qrand()%N;
	}
	return ret;
}

Mda estimate_noise_shape(DiskReadMda &X,int T,int ch) {
// Computes the expected shape of the template in a noise cluster
// which may be considered as a row in the noise covariance matrix
// X is the MxN array of raw or pre-processed data
// T is the clip size
// ch is the channel where detection takes place
	QTime timer; timer.start();
	int M=X.N1();
	int N=X.N2();
	int Tmid=(int)((T+1)/2)-1;
	int num_rand_times=10000;
	QList<double> rand_times=randsample_with_replacement(N-2*T,num_rand_times);
	for (int i=0; i<rand_times.count(); i++) rand_times[i]+=T;
	Mda rand_clips=extract_clips(X,rand_times,T);
	QList<double> peaks;
	for (int i=0; i<rand_times.count(); i++) {
		peaks << rand_clips.get(ch,Tmid,i);
	}
	Mda noise_shape; noise_shape.allocate(M,T);
	double *noise_shape_ptr=noise_shape.dataPtr();
	double *rand_clips_ptr=rand_clips.dataPtr();
	int bbb=0;
	for (int i=0; i<rand_times.count(); i++) {
		if (fabs(peaks[i])<=2) { //focus on noise clips (where amplitude is low)
			int aaa=0;
			for (int t=0; t<T; t++) {
				for (int m=0; m<M; m++) {
					noise_shape_ptr[aaa]+=peaks[i]*rand_clips_ptr[bbb];
					aaa++;
					bbb++;
				}
			}
		}
		else bbb+=M*T; //important and tricky!!!
	}
	double noise_shape_norm=0;
	for (int t=0; t<T; t++) {
		for (int m=0; m<M; m++) {
			double val=noise_shape.get(m,t);
			noise_shape_norm+=val*val;
		}
	}
	noise_shape_norm=sqrt(noise_shape_norm);
	if (noise_shape_norm) {
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				noise_shape.set(noise_shape.get(m,t)/noise_shape_norm,m,t);
			}
		}
	}
	return noise_shape;
}

Mda compute_features(Mda &clips,int num_features) {
	int M=clips.N1();
	int T=clips.N2();
	int L=clips.N3();
	Mda ret; ret.allocate(num_features,L);
	get_pca_features(M*T,L,num_features,ret.dataPtr(),clips.dataPtr());
	return ret;
}

void compute_geometric_median(int M, int N, double *output, double *input, int num_iterations=10) {
	double *weights=(double *)malloc(sizeof(double)*N);
	double *dists=(double *)malloc(sizeof(double)*N);
	for (int j=0; j<N; j++) weights[j]=1;
	for (int it=1; it<=num_iterations; it++) {
		float sumweights=0;
		for (int j=0; j<N; j++) sumweights+=weights[j];
		if (sumweights) for (int j=0; j<N; j++) weights[j]/=sumweights;
		for (int i=0; i<M; i++) output[i]=0;
		{
			//compute output
			int kkk=0;
			for (int j=0; j<N; j++) {
				int i=0;
				for (int m=0; m<M; m++) {
					output[i]+=weights[j]*input[kkk];
					i++;
					kkk++;
				}
			}
		}
		{
			//compute dists
			int kkk=0;
			for (int j=0; j<N; j++) {
				int i=0;
				double sumsqr=0;
				for (int m=0; m<M; m++) {
					double diff=output[i]-input[kkk];
					i++;
					kkk++;
					sumsqr+=diff*diff;
				}
				dists[j]=sqrt(sumsqr);
			}
		}
		{
			//compute weights
			for (int j=0; j<N; j++) {
				if (dists[j]) weights[j]=1/dists[j];
				else weights[j]=0;
			}
		}
	}
	free(dists);
	free(weights);
}

Mda compute_geometric_median_template(Mda &clips) {
	int M=clips.N1();
	int T=clips.N2();
	int L=clips.N3();
	Mda ret; ret.allocate(M,T);
	if (L==0) return ret;
	int num_features=18;
	Mda FF=compute_features(clips,num_features);
	Mda FFmm; FFmm.allocate(num_features,1);
	compute_geometric_median(num_features,L,FFmm.dataPtr(),FF.dataPtr());
	QList<double> dists;
	for (int i=0; i<L; i++) {
		double tmp=0;
		for (int a=0; a<num_features; a++) {
			double val=FF.get(a,i)-FFmm.get(a,0L);
			tmp+=val*val;
		}
		tmp=sqrt(tmp);
		dists << tmp;
	}
	QList<double> sorted_dists=dists; qSort(sorted_dists);
	double dist_cutoff=sorted_dists[(int)(sorted_dists.count()*0.3)];
	QList<int> inds;
	for (int i=0; i<dists.count(); i++) {
		if (dists[i]<=dist_cutoff) inds << i;
	}
	Mda subclips=get_subclips(clips,inds);
	return compute_mean_clip(subclips);
}

double compute_template_ip(Mda &T1,Mda &T2) {
	double ip=0;
	for (int t=0; t<T1.N2(); t++) {
		for (int m=0; m<T1.N1(); m++) {
			ip+=T1.get(m,t)*T2.get(m,t);
		}
	}
	return ip;
}

double compute_template_norm(Mda &T) {
	return sqrt(compute_template_ip(T,T));
}

QList<Subcluster> compute_subcluster_detectability_scores(DiskReadMda &pre,Mda &clips,int channel,const Remove_noise_subclusters_opts &opts) {
	QTime timer; timer.start();
	int M=clips.N1();
	int T=clips.N2();
	int L=clips.N3();
	int Tmid=(int)((T+1)/2)-1;
	QList<double> peaks;
	for (int i=0; i<L; i++) {
		peaks << clips.get(channel,Tmid,i);
	}
    Define_Shells_Opts opts2; opts2.min_shell_size=opts.min_shell_size; opts2.shell_increment=opts.shell_increment;
    QList<Shell> shells=define_shells(peaks,opts2);

    Mda noise_shape;
    #pragma omp critical
    noise_shape=estimate_noise_shape(pre,T,channel);

	QList<Subcluster> subclusters;
	for (int s=0; s<shells.count(); s++) {
		QList<int> inds_s=shells[s].inds;
		Mda clips_s=get_subclips(clips,inds_s);
		//Mda subtemplate=compute_geometric_median_template(clips_s);
		Mda subtemplate=compute_mean_clip(clips_s);
		double ip=compute_template_ip(noise_shape,subtemplate);
		Mda subtemplate_resid; subtemplate_resid.allocate(M,T);
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				subtemplate_resid.set(subtemplate.get(m,t)-ip*noise_shape.get(m,t),m,t);
			}
		}
		double subtemplate_resid_norm=compute_template_norm(subtemplate_resid);
		Subcluster SC;
		SC.inds=inds_s;
		SC.detectability_score=subtemplate_resid_norm;
		SC.peak=subtemplate.get(channel,Tmid);
		subclusters << SC;
	}
	return subclusters;
}

double compute_slope(const QList<double> &X,const QList<double> &Y) {
	if (X.count()==0) return 0;
	double mu_x=0,mu_y=0;
	for (int i=0; i<X.count(); i++) {
		mu_x+=X.value(i);
		mu_y+=Y.value(i);
	}
	mu_x/=X.count();
	mu_y/=Y.count();
	double numer=0,denom=0;
	for (int i=0; i<X.count(); i++) {
		numer+=(X.value(i)-mu_x)*(Y.value(i)-mu_y);
		denom+=(X.value(i)-mu_x)*(X.value(i)-mu_x);
	}
	if (denom) return numer/denom;
	else return 0;
}

bool remove_noise_subclusters(const QString &signal_path,const QString &firings_path,const QString &firings_out_path,const Remove_noise_subclusters_opts &opts) {

	DiskReadMda raw; raw.setPath(signal_path);
	Mda firings; firings.read(firings_path);

	int T=opts.clip_size;
	int L=firings.N2();

	QList<int> channels,labels;
	QList<double> times;
	for (int i=0; i<L; i++) {
		channels << (int)firings.get(0,i)-1; //convert to zero-based indexing
		times << firings.get(1,i) -1; //convert to zero-based indexing
		labels << (int)firings.get(2,i);
	}

	Mda clips=extract_clips(raw,times,T);

	int K=compute_max(labels);

	QList<Subcluster> subclusters;

    #pragma omp parallel for
	for (int k=1; k<=K; k++) {
		int count1=0;
		int count2=0;

        QList<int> inds_k;
        #pragma omp critical
        inds_k=find_label_inds(labels,k);

        if (inds_k.count()>0) {
            int channel;
            Mda clips_k;
            #pragma omp critical
            {
                printf("k=%d/%d\n",k,K);
                channel=channels[inds_k[0]];
                clips_k=get_subclips(clips,inds_k);
            }
			QList<Subcluster> subclusters0=compute_subcluster_detectability_scores(raw,clips_k,channel,opts);
            #pragma omp critical
            {
                for (int ii=0; ii<subclusters0.count(); ii++) {
                    Subcluster SC=subclusters0[ii];
                    QList<int> new_inds;
                    for (int j=0; j<SC.inds.count(); j++) {
                        new_inds << inds_k[SC.inds[j]];
                    }
                    SC.inds=new_inds;
                    SC.label=k;
                    subclusters << SC;
                    count1+=SC.inds.count();
                    if (SC.inds.count()>0) count2++;
                }
            }
		}
	}

	int num_subclusters=subclusters.count();
	QList<double> subcluster_abs_peaks,subcluster_scores;
	for (int ii=0; ii<num_subclusters; ii++) {
		subcluster_abs_peaks << fabs(subclusters[ii].peak);
		subcluster_scores << subclusters[ii].detectability_score;
	}
	double slope=compute_slope(subcluster_abs_peaks,subcluster_scores);
	if (slope) {
		for (int ii=0; ii<num_subclusters; ii++) {
			subcluster_scores[ii]/=slope;
			subclusters[ii].detectability_score/=slope;
		}
	}
	QList<int> labels_to_use; for (int k=0; k<=K; k++) labels_to_use << 0;
	QList<int> events_to_use; for (int i=0; i<L; i++) events_to_use << 0;
	int num_events_to_use=0;
	for (int ii=0; ii<num_subclusters; ii++) {
		if (subcluster_scores[ii]>=opts.detectability_threshold) {
			labels_to_use[subclusters[ii].label]=1;
			for (int jj=0; jj<subclusters[ii].inds.count(); jj++) {
				events_to_use[subclusters[ii].inds[jj]]=1;
				num_events_to_use++;
			}
		}
	}
	QList<int> labels_map; for (int k=0; k<=K; k++) labels_map << 0;
	int kk=1;
	for (int k=1; k<=K; k++) {
		if (labels_to_use[k]) {
			labels_map[k]=kk;
			kk++;
		}
	}
	Mda firings_out;
	firings_out.allocate(firings.N1(),num_events_to_use);
	int ee=0;
	for (int i=0; i<L; i++) {
		if (events_to_use[i]) {
			for (int jj=0; jj<firings.N1(); jj++) {
				firings_out.set(firings.get(jj,i),jj,ee);
			}
			int label=(int)firings_out.get(2,ee);
			if (label) label=labels_map[label];
			firings_out.set(label,2,ee);
			ee++;
		}
	}

	printf("Using %d/%d events in %d/%d clusters\n",num_events_to_use,L,compute_max(labels_map),K);
    firings_out.write64(firings_out_path);

	return true;
}
