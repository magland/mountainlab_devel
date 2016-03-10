#ifndef DETECT2_H
#define DETECT2_H

struct Detect_Opts {
    int clip_size;
    int sign;
    bool individual_channels;
};

bool detect2(const char *raw_path,const char *detect_path,double detect_threshold,int detect_interval,const Detect_Opts &opts);

#endif // DETECT2_H

