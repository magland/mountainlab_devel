#ifndef MVCLUSTERDETAILWIDGET_H
#define MVCLUSTERDETAILWIDGET_H

#include "diskreadmda.h"
#include "mda.h"
#include <QWidget>

class MVClusterDetailWidgetPrivate;
class MVClusterDetailWidget : public QWidget
{
	Q_OBJECT
public:
	friend class MVClusterDetailWidgetPrivate;
	MVClusterDetailWidget(QWidget *parent=0);
	virtual ~MVClusterDetailWidget();
	void setRaw(DiskReadMda &X);
	void setFirings(DiskReadMda &X);
	void setSamplingFrequency(double freq);
	void setChannelColors(const QList<QColor> &colors);
	int currentK();
	void setCurrentK(int k);
protected:
	void paintEvent(QPaintEvent *evt);
	void keyPressEvent(QKeyEvent *evt);
	void mousePressEvent(QMouseEvent *evt);
	void mouseReleaseEvent(QMouseEvent *evt);
	void mouseMoveEvent(QMouseEvent *evt);
signals:
	void signalCurrentKChanged();
private:
	MVClusterDetailWidgetPrivate *d;
};

#endif // MVCLUSTERDETAILWIDGET_H
