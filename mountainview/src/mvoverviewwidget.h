#ifndef MVOVERVIEWWIDGET_H
#define MVOVERVIEWWIDGET_H

#include <QMainWindow>
#include <QWheelEvent>
#include "diskarraymodel.h"
#include "mda.h"

class MVOverviewWidgetPrivate;

class MVOverviewWidget : public QMainWindow
{
	Q_OBJECT
public:
	friend class MVOverviewWidgetPrivate;
	explicit MVOverviewWidget(QWidget *parent = 0);
	~MVOverviewWidget();

	void setElectrodeLocations(const Mda &L);
	void setTemplates(const Mda &X);
	void setPrimaryChannels(const Mda &X);
	void setRaw(DiskArrayModel *X,bool own_it);
	void setTimesLabels(const Mda &times,const Mda &labels);
	void setCrossCorrelogramsPath(const QString &path);

	void setClips(DiskArrayModel *X,bool own_it);
	void setClipsIndex(const Mda &X);

	void updateWidgets();

protected:
	void resizeEvent(QResizeEvent *evt);

private slots:
	void slot_spike_templates_current_x_changed();
	void slot_cross_correlograms_current_unit_changed();
	void slot_statistics_widget_current_unit_changed();
	void slot_unit_activated(int num);


private:
	MVOverviewWidgetPrivate *d;

signals:

public slots:
};

#endif // MVOVERVIEWWIDGET_H
