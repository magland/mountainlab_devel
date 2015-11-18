#ifndef BANDPASS_FILTER_H
#define BANDPASS_FILTER_H

bool bandpass_filter(const char *input_path,const char *output_path,double samplefreq,double freq_min,double freq_max,double outlier_threshold);

#endif // BANDPASS_FILTER_H

