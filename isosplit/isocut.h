#ifndef isocut_h
#define isocut_h

//return true if split is statistically significant
bool isocut(int N,double &cutpoint,double *X,double threshold=1.4,int minsize=4);

#endif
