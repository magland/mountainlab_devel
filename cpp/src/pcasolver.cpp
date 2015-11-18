#include "pcasolver.h"
#include <math.h>
#include <QDebug>

class PCASolverPrivate {
public:
	int m_N1;
	int m_N2;
	float *m_vectors;
	float *m_components;
	float *m_coefficients;
	int m_component_count;
	int m_num_iterations;
	//QList<double> m_energies;
};

PCASolver::PCASolver() 
{
	d=new PCASolverPrivate;
	d->m_component_count=1;
	d->m_num_iterations=10;
	d->m_vectors=0;
	d->m_coefficients=0;
	d->m_components=0;
}
PCASolver::~PCASolver()
{
	if (d->m_coefficients) free(d->m_coefficients);
	if (d->m_components) free(d->m_components);
	delete d;
}
void PCASolver::setVectors(int N1,int N2,float *V) {
	d->m_N1=N1;
	d->m_N2=N2;
	d->m_vectors=V;
}
void PCASolver::setComponentCount(int c) {
	d->m_component_count=c;
}

float vector_norm(long N,float *v) {
	float norm=0;
	for (long j=0; j<N; j++)
		norm+=v[j]*v[j];
	norm=sqrt(norm);
	return norm;
}
void define_random_vector(long N,float *v) {
	for (long j=0; j<N; j++) {
		v[j]=(qrand()*1.0/RAND_MAX)*2-1;
	}
}
void normalize(long N,float *v) {
	float norm=vector_norm(N,v);
	if (norm==0) return;
	for (long j=0; j<N; j++)
		v[j]/=norm;
}
float inner_product(long N,float *v,float *w) {
	float ret=0;
	for (long j=0; j<N; j++)
		ret+=v[j]*w[j];
	return ret;
}
void subtract_component(long N,float *v,float *comp) {
	float ip=inner_product(N,v,comp);
	for (long j=0; j<N; j++) {
		v[j]-=ip*comp[j];
	}
}



bool PCASolver::solve() {
	int N1=d->m_N1;
	int N2=d->m_N2;
	int num_components=d->m_component_count;

	if (d->m_coefficients) free(d->m_coefficients);
	if (d->m_components) free(d->m_components);

	//Define the working vectors
	float *working_vectors=(float *)malloc(sizeof(float)*N1*N2);
	for (int vnum=0; vnum<N2; vnum++) {
		for (long j=0; j<N1; j++)
			working_vectors[vnum*N1+j]=d->m_vectors[j+N1*vnum];
	}
	
	//Define the working components

	float *working_components=(float *)malloc(sizeof(float)*N1*num_components);
	for (int cnum=0; cnum<num_components; cnum++) {
		for (long j=0; j<N1; j++)
			working_components[cnum*N1+j]=0;
	}
	d->m_components=working_components;
	
	//allocate the coefficients and energies
	float *coeffs=(float *)malloc(sizeof(float)*N2*num_components);
	d->m_coefficients=coeffs;
	//double *energies=(double *)malloc(sizeof(double)*num_components);
	
	for (int current_component=0; current_component<num_components; current_component++) {
		float *component_vector=&working_components[current_component*N1];
		float component_norm=vector_norm(N1,component_vector);
		if (component_norm<0.1) {
			define_random_vector(N1,component_vector);
			normalize(N1,component_vector);
		}
		for (long it=0; it<d->m_num_iterations; it++) {	
			float *hold=(float *)malloc(sizeof(float)*N1);
			if (!hold) {
				qWarning() << "Unable to allocate hold of size:" << N1;
			}
			for (int j=0; j<N1; j++) hold[j]=0;
			for (int vnum=0; vnum<N2; vnum++) {
				float ip=0;
				ip=inner_product(N1,&working_vectors[vnum*N1],component_vector);
				for (int j=0; j<N1; j++) {
					hold[j]+=ip*working_vectors[vnum*N1+j];
				}
			}
			normalize(N1,hold);
			for (int j=0; j<N1; j++) component_vector[j]=hold[j];
			free(hold);
		}
		//Compute coefficients
		for (long vnum=0; vnum<N2; vnum++) {
			float ip0=0;
			ip0=inner_product(N1,&working_vectors[vnum*N1],&working_components[current_component*N1]);
			coeffs[vnum+N2*current_component]=ip0;
		}
		//Compute energy (lambda)
		double val=0;
		for (int vnum=0; vnum<N2; vnum++) {
			val+=coeffs[vnum+N2*current_component]*coeffs[vnum+N2*current_component];
		}
		//energies[current_component]=val;
		//Subtract this component from the working vectors
		for (int vnum=0; vnum<N2; vnum++) {
			subtract_component(N1,&working_vectors[N1*vnum],&working_components[N1*current_component]);
		}
	}
	

	free(working_vectors);
	
	return true;
}
float *PCASolver::components() {
	return d->m_components;
}
float *PCASolver::coefficients() {
	return d->m_coefficients;
}
void PCASolver::setNumIterations(int val) {
	d->m_num_iterations=val;
}
//QList<double> PCASolver::energies() const {
//	return d->m_energies;
//}
//void PCASolver::setWeights(const QList<float> &W) {
//	d->m_weights.allocate(W.count(),1);
//	d->m_weights.setDataX(W,0);
//}
