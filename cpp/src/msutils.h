#ifndef MSUTILS_H
#define MSUTILS_H

#include <QList>
#include "diskreadmda.h"

double compute_min(const QList<double> &X);
double compute_max(const QList<double> &X);
int compute_max(const QList<int> &X);
Mda extract_clips(DiskReadMda &X,const QList<long> &times,int clip_size);
Mda extract_clips(DiskReadMda &X,const QList<long> &times,const QList<int> &channels,int clip_size);
Mda compute_mean_clip(Mda &clips);

#endif // MSUTILS_H

