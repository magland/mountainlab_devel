#ifndef MVSTATISTICSWIDGET_H
#define MVSTATISTICSWIDGET_H

#include <QWidget>
#include "mda.h"
#include "diskarraymodel.h"

class MVStatisticsWidgetPrivate;
class MVStatisticsWidget : public QWidget
{
	Q_OBJECT
public:
	friend class MVStatisticsWidgetPrivate;
	MVStatisticsWidget();
	virtual ~MVStatisticsWidget();

	void setTimesLabels(const Mda &times,const Mda &labels);
	void setPrimaryChannels(const Mda &primary_channels);
	void setRaw(DiskArrayModel *X);
	void updateStatistics();

	int currentUnit();
	void setCurrentUnit(int unit);

signals:
	void currentUnitChanged();

private slots:
	void slot_item_clicked();

private:
	MVStatisticsWidgetPrivate *d;
};

#endif // MVSTATISTICSWIDGET_H
