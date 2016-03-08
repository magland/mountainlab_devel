#ifndef MVCLUSTERVIEW_H
#define MVCLUSTERVIEW_H

#include <QWidget>
#include "mda.h"
#include "mvutils.h"
#include "affinetransformation.h"

#define MVCV_MODE_HEAT_DENSITY 1
#define MVCV_MODE_LABEL_COLORS 2

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
	void setMode(int mode);
	void setCurrentEvent(MVEvent evt,bool do_emit=false);
	MVEvent currentEvent();
	AffineTransformation transformation();
	void setTransformation(const AffineTransformation &T);
signals:
	void currentEventChanged();
	void transformationChanged();
protected:
	void paintEvent(QPaintEvent *evt);
	void mouseMoveEvent(QMouseEvent *evt);
	void mousePressEvent(QMouseEvent *evt);
	void mouseReleaseEvent(QMouseEvent *evt);
	void wheelEvent(QWheelEvent *evt);
private slots:
	void slot_emit_transformation_changed();
private:
	MVClusterViewPrivate *d;
};

#endif // MVCLUSTERVIEW_H
