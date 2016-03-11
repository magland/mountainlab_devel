#ifndef MSUTILS_H
#define MSUTILS_H

#include <QList>
#include "diskreadmda.h"

double compute_min(const QList<double> &X);
double compute_max(const QList<double> &X);
double compute_max(long N,double *X);
int compute_max(const QList<int> &X);
Mda extract_clips(DiskReadMda &X,const QList<double> &times,int clip_size);
Mda extract_clips(DiskReadMda &X,const QList<long> &times,int clip_size);
Mda extract_clips(DiskReadMda &X,const QList<long> &times,const QList<int> &channels,int clip_size);
Mda compute_mean_clip(Mda &clips);
double compute_mean(const QList<double> &X);
double compute_stdev(const QList<double> &X);
Mda grab_clips_subset(Mda &clips,const QList<int> &inds);

#endif // MSUTILS_H

