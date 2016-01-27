#ifndef ISOBRANCH_H
#define ISOBRANCH_H

bool isobranch(const char *clips_input_path,const char *labels_output_path,int min_cluster_size=500,int num_features=3,float isocut_threshold=1.2,int K_init=30);

#endif // ISOBRANCH_H

