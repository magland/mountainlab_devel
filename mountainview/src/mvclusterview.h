#ifndef MVCLUSTERVIEW_H
#define MVCLUSTERVIEW_H

#include <QWidget>
#include "mda.h"

class MVClusterViewPrivate;
class MVClusterView : public QWidget
{
	Q_OBJECT
public:
	friend class MVClusterViewPrivate;
	MVClusterView(QWidget *parent=0);
	virtual ~MVClusterView();
	void setData(const Mda &X);
protected:
	void paintEvent(QPaintEvent *evt);
private:
	MVClusterViewPrivate *d;
};

#endif // MVCLUSTERVIEW_H
