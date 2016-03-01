#ifndef REMOVE_NOISE_SUBCLUSTERS_H
#define REMOVE_NOISE_SUBCLUSTERS_H

struct Remove_noise_subclusters_opts {
	int clip_size;
	double detectability_threshold;
	double shell_increment;
	double min_shell_size;
};

bool remove_noise_subclusters(const char *pre_path,const char *firings_path,const char *firings_out_path,const Remove_noise_subclusters_opts &opts);

#endif // REMOVE_NOISE_SUBCLUSTERS_H

