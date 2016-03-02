#ifndef MVCLUSTERDETAILWIDGET_H
#define MVCLUSTERDETAILWIDGET_H

#include "diskreadmda.h"
#include "mda.h"
#include <QWidget>
#include <QScrollArea>

class MVClusterDetailWidgetPrivate;
class MVClusterDetailWidget : public QWidget
{
	Q_OBJECT
public:
	friend class MVClusterDetailWidgetPrivate;
	friend class MVClusterDetailWidgetScrollArea;
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
	void wheelEvent(QWheelEvent *evt);
signals:
	void signalCurrentKChanged();
	void signalZoomedIn();
private:
	MVClusterDetailWidgetPrivate *d;
};

class MVClusterDetailWidgetScrollArea : public QScrollArea {
	Q_OBJECT
public:
	MVClusterDetailWidgetScrollArea(QWidget *parent=0);
	void setTheWidget(MVClusterDetailWidget *W);
	MVClusterDetailWidget *theWidget();
private slots:
	void slot_current_k_changed();
	void slot_zoomed_in();
private:
	MVClusterDetailWidget *the_widget;
	void ensure_visible(double x);
};

#endif // MVCLUSTERDETAILWIDGET_H
