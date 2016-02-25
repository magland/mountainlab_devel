#ifndef isosplit2_h
#define isosplit2_h

#include "mda.h"
#include <QList>

QList<int> isosplit2(Mda &X,float isocut_threshold=1.5,int K_init=30);

#endif
