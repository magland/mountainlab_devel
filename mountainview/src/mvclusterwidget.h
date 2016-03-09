#ifndef MVCLUSTERWIDGET_H
#define MVCLUSTERWIDGET_H

#include <QWidget>
#include "mvutils.h"
#include "affinetransformation.h"

class MVClusterWidgetPrivate;
class MVClusterWidget : public QWidget
{
	Q_OBJECT
public:
	friend class MVClusterWidgetPrivate;
	MVClusterWidget();
	virtual ~MVClusterWidget();
	void setData(const Mda &X);
	void setTimes(const QList<double> &times);
	void setLabels(const QList<int> &labels);
	void setOutlierScores(const QList<double> &outlier_scores);
	void setCurrentEvent(const MVEvent &evt);
	void setClipSize(int clip_size);
	void setRaw(const DiskReadMda &X);
    void setTransformation(const AffineTransformation &T);
    void setDensityViewVisible(bool val);
    void setColorViewVisible(bool val);
    void setClipsViewVisible(bool val);
	MVEvent currentEvent();

signals:
	void currentEventChanged();
private slots:
	void slot_view_current_event_changed();
	void slot_view_transformation_changed();
private:
	MVClusterWidgetPrivate *d;
};

#endif // MVCLUSTERWIDGET_H

