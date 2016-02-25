#ifndef BRANCH_CLUSTER_V1_H
#define BRANCH_CLUSTER_V1_H

struct Branch_Cluster_Opts {
    int clip_size;
    int min_section_count;
    double section_increment;
    int num_features;
};

bool branch_cluster(const char *raw_path,const char *detect_path,const char *adjacency_matrix_path,const char *output_firings_path,const Branch_Cluster_Opts &opts);

#endif // BRANCH_CLUSTER_V1_H

