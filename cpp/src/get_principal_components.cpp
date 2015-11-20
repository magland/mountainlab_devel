#include "get_principal_components.h"
//Depends on libpca which depends on armadillo
#include "pca.h"
#include <QDebug>

#include "pcasolver.h"
void get_principal_components_2(int M,int N,int ncomp,float *components,float *data) {
    PCASolver PCA;
    PCA.setVectors(M,N,data);
    PCA.setComponentCount(ncomp);
    PCA.setNumIterations(100);
    PCA.solve();
    float *components0=PCA.components();
    for (int i=0; i<ncomp*M; i++) {
        components[i]=components0[i];
    }
}

void get_principal_components(int M,int N,int ncomp,float *components,float *data) {
    const int nvar = M;
    stats::pca pca(nvar);

    for (int i=0; i<N; i++) {
        std::vector<double> V;
        for (int m=0; m<M; m++) {
            V.push_back(data[m+M*i]);
        }
        pca.add_record(V);
    }

    pca.solve();

    for (int cc=0; cc<ncomp; cc++) {
        std::vector<double> prinvec = pca.get_principal(cc);
        //for (int m=0; m<M; m++) {
        //    components[m+cc*M]=prinvec[m];
        //}
        for (int n=0; n<N; n++) {
            float val=prinvec[n];
            components[cc+n*ncomp]=val;
        }
    }
}

