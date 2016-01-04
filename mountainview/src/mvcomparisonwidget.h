#ifndef MVCOMPARISONWIDGET_H
#define MVCOMPARISONWIDGET_H

#include <QWidget>
#include <QThread>
#include <QWheelEvent>
#include "diskarraymodel.h"
#include "mda.h"

class MVComparisonWidgetPrivate;

class MVComparisonWidget : public QWidget
{
	Q_OBJECT
public:
	friend class MVComparisonWidgetPrivate;
	explicit MVComparisonWidget(QWidget *parent = 0);
	~MVComparisonWidget();

	void setElectrodeLocations(const Mda &L);
	void setTemplates(const Mda &X);
	void setPrimaryChannels(const Mda &X);
	void setRaw(DiskArrayModel *X,bool own_it);
	void setTimesLabels(const Mda &times,const Mda &labels);
	void setCrossCorrelogramsPath(const QString &path);

	void setClips(const QList<DiskArrayModel *> &C,bool own_it);
	void setUnitNumbers(const QList<int> &numbers);

	void updateWidgets();

signals:
	void currentClipNumberChanged();

protected:
	void resizeEvent(QResizeEvent *evt);

private slots:
	void slot_compute_templates();
	void slot_clips_view_current_x_changed();
	void slot_selected_data_points_changed();


private:
	MVComparisonWidgetPrivate *d;

signals:

public slots:
};



#endif // MVCOMPARISONWIDGET_H
