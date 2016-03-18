#ifndef BRANCH_CLUSTER_V2_H
#define BRANCH_CLUSTER_V2_H

struct Branch_Cluster_V2_Opts {
    int clip_size;
	int min_shell_size;
	double shell_increment;
    int num_features;
	int detect_interval;
};

bool branch_cluster_v2(const char *raw_path,const char *detect_path,const char *adjacency_matrix_path,const char *output_firings_path,const Branch_Cluster_V2_Opts &opts);

#endif // BRANCH_CLUSTER_V2_H

