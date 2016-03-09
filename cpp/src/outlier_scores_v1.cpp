#include "outlier_scores_v1.h"

#include "mda.h"
#include "diskreadmda.h"
#include "msutils.h"
#include "remove_noise_subclusters.h"
#include <math.h>

QList<double> grab_sublist(const QList<double> &X,const QList<int> &inds);
QList<double> compute_outlier_scores(Mda &clips,Mda &random_clips);

bool outlier_scores_v1(const char *raw_path, const char *firings_in_path, const char *firings_out_path, const Outlier_Scores_Opts &opts)
{
    Mda firings; firings.read(firings_in_path);
    DiskReadMda X; X.setPath(raw_path);
    int N=X.N2();
    int L=firings.N2();
    QList<double> times;
    QList<int> labels;
    QList<double> peaks;
    QList<double> scores;
    for (int i=0; i<L; i++) {
        times << firings.value(1,i);
        labels << (int)firings.value(2,i);
        peaks << firings.value(3,i);
        scores << 0;
    }

    int interval=(N/5000);
    QList<double> ttt; for (int i=0; i<N; i+=interval) ttt << i;
    Mda random_clips=extract_clips(X,ttt,opts.clip_size);

    Define_Shells_Opts opts2; opts2.min_shell_size=opts.min_shell_size; opts2.shell_increment=opts.shell_increment;

    int K=compute_max(labels);
    for (int k=1; k<=K; k++) {
        QList<int> inds_k;
        for (int i=0; i<L; i++) {
            if (labels[i]==k) inds_k << i;
        }
        QList<double> peaks_k=grab_sublist(peaks,inds_k);

        QList<Shell> shells=define_shells(peaks_k,opts2);
        for (int s=0; s<shells.count(); s++) {
            QList<int> inds_ks;
            for (int j=0; j<shells[s].inds.count(); j++) {
                inds_ks << inds_k[shells[s].inds[j]];
            }
            QList<double> times_ks=grab_sublist(times,inds_ks);
            Mda clips_ks=extract_clips(X,times_ks,opts.clip_size);
            QList<double> scores_ks=compute_outlier_scores(clips_ks,random_clips);
            for (int j=0; j<inds_ks.count(); j++) {
                scores[inds_ks[j]]=scores_ks[j];
            }
        }
    }

    int R=firings.N1();
    if (R<5) R=5;
    Mda firings2; firings2.allocate(R,L);
    for (int i=0; i<L; i++) {
        for (int r=0; r<R; r++) {
            firings2.setValue(firings.value(r,i),r,i);
        }
        firings2.setValue(scores[i],4,i);
    }
    firings2.write(firings_out_path);

    return true;
}
Mda get_template_weights(Mda &template0,int num_pix) {
    int M=template0.N1();
    int T=template0.N2();
    Mda ret; ret.allocate(M,T);
    for (int m=0; m<M; m++) {
        for (int t=0; t<T; t++) {
            double val=0;
            for (int dt=-T; dt<=T; dt++) {
                if ((0<=t+dt)&&(t+dt<T)) {
                    val+=fabs(template0.value(m,t+dt))*exp(-0.5*dt*dt/(num_pix*num_pix));
                }
            }
            ret.setValue(val,m,t);
        }
    }
    return ret;
}

QList<double> compute_outlier_scores(Mda &clips,Mda &random_clips) {
    int M=clips.N1();
    int T=clips.N2();
    int L=clips.N3();
    int num_random_clips=random_clips.N3();
    Mda template0=compute_mean_clip(clips);
    Mda weights=get_template_weights(template0,6);
    Mda random_clips_weighted; random_clips_weighted.allocate(M,T,num_random_clips);
    {
        int aaa=0;
        for (int i=0; i<num_random_clips; i++) {
            int bbb=0;
            for (int t=0; t<T; t++) {
                for (int m=0; m<M; m++) {
                    random_clips_weighted.setValue1(random_clips.value1(aaa)*weights.value1(bbb),aaa);
                    aaa++; bbb++;
                }
            }
        }
    }
    Mda clips_weighted; clips_weighted.allocate(M,T,L);
    {
        int aaa=0;
        for (int i=0; i<L; i++) {
            int bbb=0;
            for (int t=0; t<T; t++) {
                for (int m=0; m<M; m++) {
                    clips_weighted.setValue1(clips.value1(aaa)*weights.value1(bbb),aaa);
                    aaa++; bbb++;
                }
            }
        }
    }
    Mda template_weighted=template0;
    for (int t=0; t<T; t++) {
        for (int m=0; m<M; m++) {
            template_weighted.setValue(template0.value(m,t)*weights.value(m,t),m,t);
        }
    }
    Mda diffs_weighted; diffs_weighted.allocate(M,T,L);
    {
        int aaa=0;
        for (int i=0; i<L; i++) {
            int bbb=0;
            for (int t=0; t<T; t++) {
                for (int m=0; m<M; m++) {
                    diffs_weighted.setValue1(clips_weighted.value1(aaa)-template_weighted.value1(bbb),aaa);
                    aaa++; bbb++;
                }
            }
        }
    }
    QList<double> vals1;
    {
        int aaa=0;
        for (int i=0; i<L; i++) {
            double val=0;
            for (int t=0; t<T; t++) {
                for (int m=0; m<M; m++) {
                    val+=diffs_weighted.value1(aaa)*diffs_weighted.value1(aaa);
                    aaa++;
                }
            }
            vals1 << val;
        }
    }
    QList<double> vals2;
    {
        int aaa=0;
        for (int i=0; i<num_random_clips; i++) {
            double val=0;
            for (int t=0; t<T; t++) {
                for (int m=0; m<M; m++) {
                    val+=random_clips_weighted.value1(aaa)*random_clips_weighted.value1(aaa);
                    aaa++;
                }
            }
            vals2 << val;
        }
    }
    double mu0=compute_mean(vals2);
    double sigma0=compute_stdev(vals2);
    if (!sigma0) sigma0=1;
    QList<double> scores;
    for (int i=0; i<L; i++) {
        scores << (vals1[i]-mu0)/sigma0;
    }

    return scores;
}

QList<double> grab_sublist(const QList<double> &X,const QList<int> &inds) {
    QList<double> ret;
    for (int i=0; i<inds.count(); i++) ret << X[inds[i]];
    return ret;
}
