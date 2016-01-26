#ifndef isosplit_h
#define isosplit_h

#include "mda.h"
#include <QVector>

QVector<int> isosplit(Mda &X,float ks_threshold,int K_init);

#endif
