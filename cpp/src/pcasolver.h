#ifndef pcasolver_H
#define pcasolver_H

class PCASolverPrivate;
class PCASolver {
public:
	friend class PCASolverPrivate;
	PCASolver();
	virtual ~PCASolver();
	void setVectors(int N1,int N2,float *V);
	void setComponentCount(int c);
	float *components();
	float *coefficients();
	void setNumIterations(int val=10);
	///float *energies();
	bool solve();
private:
	PCASolverPrivate *d;
};

float vector_norm(long N,float *v);
void define_random_vector(long N,float *v);
void normalize(long N,float *v);
float inner_product(long N,float *v,float *w);
void subtract_component(long N,float *v,float *comp);

#endif
