#include "detect2.h"
#include "diskreadmda.h"
#include <math.h>

QList<double> do_detect2(const QList<double> &vals,int detect_interval,double detect_threshold,int clip_size);

bool detect2(const char *raw_path, const char *detect_path, double detect_threshold, int detect_interval, const Detect_Opts &opts)
{
    DiskReadMda X; X.setPath(raw_path);
    int M=X.N1();
    int N=X.N2();

    QList<int> channels;
    QList<double> timepoints;
    if (opts.individual_channels) {
        for (int m=0; m<M; m++) {
            QList<double> vals;
            for (int n=0; n<N; n++) vals << X.value(m,n);
            if (opts.sign>0) {
                for (int n=0; n<N; n++) if (vals[n]<0) vals[n]=0;
            }
            if (opts.sign<0) {
                for (int n=0; n<N; n++) if (vals[n]>0) vals[n]=0;
            }
            QList<double> times=do_detect2(vals,detect_interval,detect_threshold,opts.clip_size);
            for (int a=0; a<times.count(); a++) {
                timepoints << times[a]+1;
                channels << m+1;
            }
        }
    }
    else {
        QList<double> vals;
        for (int n=0; n<N; n++) {
            int best_m=0;
            double maxabsval=fabs(X.value(0,n));
            for (int m=0; m<M; m++) {
                double val=X.value(m,n);
                if (fabs(val)>maxabsval) {
                    maxabsval=fabs(val);
                    best_m=m;
                }
            }
            vals << X.value(best_m,n);
        }
        if (opts.sign>0) {
            for (int n=0; n<N; n++) if (vals[n]<0) vals[n]=0;
        }
        if (opts.sign<0) {
            for (int n=0; n<N; n++) if (vals[n]>0) vals[n]=0;
        }
        QList<double> times=do_detect2(vals,detect_interval,detect_threshold,opts.clip_size);
        for (int a=0; a<times.count(); a++) {
            timepoints << times[a]+1;
            channels << 0;
        }
    }

    int L=timepoints.count();
    Mda detect; detect.allocate(2,L);
    for (int i=0; i<L; i++) {
        detect.setValue(channels[i],0,i);
        detect.setValue(timepoints[i],1,i);
    }
    detect.write64(detect_path);

    return true;
}

QList<double> do_detect2(const QList<double> &vals,int detect_interval,double detect_threshold,int clip_size) {
    int N=vals.count();
    int *to_use=(int *)malloc(sizeof(int)*N);
    for (int n=0; n<N; n++) to_use[n]=0;
    int last_best_ind=0;
    double last_best_abs_val=0;
    for (int n=clip_size; n<N-clip_size; n++) {
        double val=vals[n];
        if (n-last_best_ind>detect_interval) last_best_abs_val=0;
        if (fabs(val)>=detect_threshold) {
            if (last_best_abs_val>0) {
                if (fabs(val)>last_best_abs_val) {
                    to_use[n]=1;
                    to_use[last_best_ind]=0;
                    last_best_ind=n;
                    last_best_abs_val=fabs(val);
                }
            }
            else {
                to_use[n]=1;
                last_best_ind=n;
                last_best_abs_val=fabs(val);
            }
        }
    }
    QList<double> times;
    for (int n=0; n<N; n++) {
        if (to_use[n]) {
            times << n;
        }
    }
    return times;
}
