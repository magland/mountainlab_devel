#ifndef OUTLIER_SCORES_V1_H
#define OUTLIER_SCORES_V1_H

struct Outlier_Scores_Opts {
    int clip_size;
    int min_shell_size;
    double shell_increment;
};

bool outlier_scores_v1(const char *raw_path,const char *firings_in_path,const char *firings_out_path,const Outlier_Scores_Opts &opts);

#endif // OUTLIER_SCORES_V1_H

