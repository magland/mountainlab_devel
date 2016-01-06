#ifndef CVWIDGET_H
#define CVWIDGET_H

#include <QWidget>
#include "diskreadmda.h"

class CVWidgetPrivate;
class CVWidget : public QWidget
{
	Q_OBJECT
public:
	friend class CVWidgetPrivate;
	explicit CVWidget(QWidget *parent = 0);
	~CVWidget();

	Q_INVOKABLE void setFeatures(const DiskReadMda &F);
	Q_INVOKABLE void setClips(const DiskReadMda &X);
	Q_INVOKABLE void setLabels(const DiskReadMda &L);
	Q_INVOKABLE void setRange(float xmin,float xmax,float ymin,float ymax,float zmin,float zmax);
	Q_INVOKABLE void autoSetRange();
	Q_INVOKABLE void refresh();

	void setSelectedDataPointIndices(const QList<int> &L);
	QList<int> selectedDataPointIndices();
	void setNumDataPointsToSelect(int num);
    void setLabelStrings(const QStringList &strings);

signals:
	void selectedDataPointsChanged();
private:
	CVWidgetPrivate *d;

signals:

public slots:
};

#endif // CVWIDGET_H
