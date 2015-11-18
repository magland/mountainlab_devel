#include "get_principal_components.h"
#include "pcasolver.h"

void get_principal_components(int M,int N,int ncomp,float *components,float *data) {
    PCASolver PCA;
    PCA.setVectors(M,N,data);
    PCA.setComponentCount(ncomp);
    PCA.solve();
    float *components0=PCA.components();
    for (int i=0; i<ncomp*M; i++) {
        components[i]=components0[i];
    }
}

