#ifndef MVUTILS_H
#define MVUTILS_H

#include "diskarraymodel.h"

Mda compute_mean_waveform(DiskArrayModel *C);
Mda compute_features(DiskArrayModel *C);
Mda compute_features(const QList<DiskArrayModel *> &C);

#endif // MVUTILS_H

