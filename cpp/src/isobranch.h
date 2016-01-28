#ifndef ISOBRANCH_H
#define ISOBRANCH_H

#include <QList>

bool isobranch(
        const char *clips_input_path,
        const char *labels_output_path,
        const QList<float> &branch_thresholds,
        int num_features=3,
        float isocut_threshold=1.2,
        int K_init=30
);

#endif // ISOBRANCH_H

