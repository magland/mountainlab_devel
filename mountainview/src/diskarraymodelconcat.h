#ifndef DISKARRAYMODELCONCAT_H
#define DISKARRAYMODELCONCAT_H

#include "diskarraymodel.h"

class DiskArrayModelConcatPrivate;
class DiskArrayModelConcat : public DiskArrayModel
{
	Q_OBJECT
public:
	friend class DiskArrayModelConcatPrivate;
	DiskArrayModelConcat();
	virtual ~DiskArrayModelConcat();

	virtual Mda loadData(int scale,int t1,int t2);
	virtual int size(int dim);
	virtual int dim3();

	void setRange(int t1,int t2);
signals:

private slots:

private:
	DiskArrayModelConcatPrivate *d;
};

#endif // DISKARRAYMODELCONCAT_H
