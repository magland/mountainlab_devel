#ifndef MVCLUSTERVIEW_H
#define MVCLUSTERVIEW_H

#include <QWidget>
#include "mda.h"
#include "mvutils.h"

class MVClusterViewPrivate;
class MVClusterView : public QWidget
{
	Q_OBJECT
public:
	friend class MVClusterViewPrivate;
	MVClusterView(QWidget *parent=0);
	virtual ~MVClusterView();
	void setData(const Mda &X);
	void setTimes(const QList<double> &times);
	void setLabels(const QList<int> &labels);
	void setCurrentEvent(MVEvent evt);
	MVEvent currentEvent();
signals:
	void currentEventChanged();
protected:
	void paintEvent(QPaintEvent *evt);
	void mouseMoveEvent(QMouseEvent *evt);
	void mousePressEvent(QMouseEvent *evt);
	void mouseReleaseEvent(QMouseEvent *evt);
	void wheelEvent(QWheelEvent *evt);

private:
	MVClusterViewPrivate *d;
};

#endif // MVCLUSTERVIEW_H
