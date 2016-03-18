#include "get_principal_components.h"
#include <math.h>
#include <QDebug>
#ifdef QT_VERSION
#include <QApplication>
#include <QTime>
#endif

#ifdef QT_VERSION
static QTime s_timer;
#endif
void process_gui_events_if_needed() {
#ifdef QT_VERSION
	if (s_timer.elapsed()>500) {
		qApp->processEvents();
		s_timer.restart();
	}
#endif
}

void make_random_vector(int M,double *v) {
	for (int i=0; i<M; i++) v[i]=(qrand()*1.0/RAND_MAX)*2-1;
}

void XXT_vector_mult(int M,int N,double *A,double *v) {
	double *v2=(double *)malloc(sizeof(double)*N);
	for (int n=0; n<N; n++) v2[n]=0;
	int ii=0;
	for (int n=0; n<N; n++) {
		for (int m=0; m<M; m++) {
			v2[n]+=A[ii]*v[m];
			ii++;
		}
	}
	ii=0;
	for (int m=0; m<M; m++) v[m]=0;
	for (int n=0; n<N; n++) {
		for (int m=0; m<M; m++) {
			v[m]+=A[ii]*v2[n];
			ii++;
		}
	}
	free(v2);
}

void normalize_vector(int M,double *v) {
	double sumsqr=0;
	for (int m=0; m<M; m++) sumsqr+=v[m]*v[m];
	if (sumsqr>0) {
		double norm=sqrt(sumsqr);
		for (int m=0; m<M; m++) v[m]/=norm;
	}
}

void get_top_component(int M,int N,double *comp,double *data,int num_iterations=30) {
	make_random_vector(M,comp);
	for (int ii=0; ii<num_iterations; ii++) {
		process_gui_events_if_needed();
		XXT_vector_mult(M,N,data,comp);
		normalize_vector(M,comp);
	}
}

void subtract_component_from_vector(int M,double *v,double *data) {
	double ip=0;
	for (int m=0; m<M; m++) ip+=v[m]*data[m];
	for (int m=0; m<M; m++) data[m]-=ip*v[m];
}

void subtract_component_from_data(int M,int N,double *v,double *data) {
	for (int n=0; n<N; n++) {
		subtract_component_from_vector(M,v,&data[M*n]);
	}
}

double compute_inner_product(int M,double *v1,double *v2) {
	double ip=0;
	for (int m=0; m<M; m++) ip+=v1[m]*v2[m];
	return ip;
}

void do_pca(int M,int N,int ncomp,double *components,double *features,double *data) {
	#ifdef QT_VERSION
	s_timer.start();
	#endif
	double *working_data=(double *)malloc(sizeof(double)*M*N);
	for (int ii=0; ii<M*N; ii++) working_data[ii]=data[ii];
	double *v=(double *)malloc(sizeof(double)*M);
	for (int k=0; k<ncomp; k++) {
		process_gui_events_if_needed();
		get_top_component(M,N,v,working_data);
		subtract_component_from_data(M,N,v,working_data);
		memcpy(&components[M*k],v,M*sizeof(double));
	}
	for (int n=0; n<N; n++) {
		for (int k=0; k<ncomp; k++) {
			features[k+ncomp*n]=compute_inner_product(M,&components[M*k],&data[M*n]);
		}
	}
}

/*
void get_pca_features_old(int M,int N,int ncomp,double *features,double *data) {
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
		for (int n=0; n<N; n++) {
			float val=prinvec[n];
			features[cc+n*ncomp]=val;
		}
	}
}
*/


void get_pca_components(int M, int N, int ncomp, float *components, float *data)
{
	double *data2=(double *)malloc(sizeof(double)*M*N);
	double *features2=(double *)malloc(sizeof(double)*N*ncomp);
	double *components2=(double *)malloc(sizeof(double)*M*ncomp);
	for (int ii=0; ii<M*N; ii++) data2[ii]=data[ii];
	do_pca(M,N,ncomp,components2,features2,data2);
	for (int ii=0; ii<M*ncomp; ii++) components[ii]=components2[ii];
	free(data2);
	free(features2);
	free(components2);
}
void get_pca_features(int M,int N,int ncomp,float *features,float *data) {
	double *data2=(double *)malloc(sizeof(double)*M*N);
	double *features2=(double *)malloc(sizeof(double)*N*ncomp);
	double *components2=(double *)malloc(sizeof(double)*M*ncomp);
	for (int ii=0; ii<M*N; ii++) data2[ii]=data[ii];
	do_pca(M,N,ncomp,components2,features2,data2);
	for (int ii=0; ii<N*ncomp; ii++) features[ii]=features2[ii];
	free(data2);
	free(features2);
	free(components2);
}

void get_pca_features(int M, int N, int ncomp, double *features, double *data)
{
	double *components2=(double *)malloc(sizeof(double)*M*ncomp);
	do_pca(M,N,ncomp,components2,features,data);
	free(components2);
}


/*




//Depends on libpca which depends on armadillo
#include "pca.h"

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

void get_pca_features_old(int M,int N,int ncomp,float *features,float *data) {
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
        for (int n=0; n<N; n++) {
            float val=prinvec[n];
			features[cc+n*ncomp]=val;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

*/