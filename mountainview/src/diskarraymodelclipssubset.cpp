#include "diskarraymodelclipssubset.h"
#include <QDebug>

class DiskArrayModelClipsSubsetPrivate {
public:
	DiskArrayModelClipsSubset *q;
	int m_t1,m_t2;
};

DiskArrayModelClipsSubset::DiskArrayModelClipsSubset() : DiskArrayModel()
{
	d=new DiskArrayModelClipsSubsetPrivate;
	d->q=this;
	d->m_t1=0;
	d->m_t2=0;
}

DiskArrayModelClipsSubset::~DiskArrayModelClipsSubset()
{
	delete d;
}

Mda DiskArrayModelClipsSubset::loadData(int scale, int t1, int t2)
{
	int t1_b=(t1*scale+d->m_t1)/scale;
	int t2_b=(t2*scale+d->m_t1)/scale;
	return DiskArrayModel::loadData(scale,t1_b,t2_b);
}

int DiskArrayModelClipsSubset::size(int dim)
{
	if (dim==0) return DiskArrayModel::size(dim);
	else if (dim==1) return d->m_t2-d->m_t1;
	else return 1;
}

int DiskArrayModelClipsSubset::dim3()
{
	int dim3a=DiskArrayModel::dim3();
	int size1a=DiskArrayModel::size(1);
	int size1=size(1);

	if (size1==0) return dim3a;
	if (size1a==0) return dim3a;
	return (int)(dim3a*(size1*1.0/size1a));
}

void DiskArrayModelClipsSubset::setRange(int t1, int t2)
{
	d->m_t1=t1;
	d->m_t2=t2;
}