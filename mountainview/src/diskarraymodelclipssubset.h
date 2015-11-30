#ifndef DISKARRAYMODELCLIPSSUBSET_H
#define DISKARRAYMODELCLIPSSUBSET_H

#include "diskarraymodel.h"

class DiskArrayModelClipsSubsetPrivate;
class DiskArrayModelClipsSubset : public DiskArrayModel
{
	Q_OBJECT
public:
	friend class DiskArrayModelClipsSubsetPrivate;
	DiskArrayModelClipsSubset();
	virtual ~DiskArrayModelClipsSubset();

	virtual Mda loadData(int scale,int t1,int t2);
	virtual int size(int dim);
	virtual int dim3();

	void setRange(int t1,int t2);
signals:

private slots:

private:
	DiskArrayModelClipsSubsetPrivate *d;
};

#endif // DISKARRAYMODELCLIPSSUBSET_H
