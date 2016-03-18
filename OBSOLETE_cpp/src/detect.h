#ifndef DETECT_H
#define DETECT_H

bool detect(const char *input_path,const char *output_path,int inner_window_width,int outer_window_width,double threshold,bool normalize,bool individual_channels,int clip_size);

#endif // DETECT_H

