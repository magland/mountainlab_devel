#ifndef MVUNITWIDGET_H
#define MVUNITWIDGET_H

#include <QMainWindow>
#include <QWheelEvent>
#include "diskarraymodel.h"
#include "mda.h"

class MVUnitWidgetPrivate;

class MVUnitWidget : public QMainWindow
{
	Q_OBJECT
public:
	friend class MVUnitWidgetPrivate;
	explicit MVUnitWidget(QWidget *parent = 0);
	~MVUnitWidget();

	void setElectrodeLocations(const Mda &L);
	void setTemplates(const Mda &X);
	void setTemplatesWhitened(const Mda &X);
	void setPrimaryChannels(const Mda &X);
	void setRaw(DiskArrayModel *X,bool own_it);
	void setRawWhitened(DiskArrayModel *X,bool own_it);
	void setTimesLabels(const Mda &times,const Mda &labels);
	void setCrossCorrelogramsPath(const QString &path);

	void setClips(DiskArrayModel *C,bool own_it);
	void setUnitNumber(int num);

	void updateWidgets();

protected:
	void resizeEvent(QResizeEvent *evt);

private slots:


private:
	MVUnitWidgetPrivate *d;

signals:

public slots:
};

#endif // MVUNITWIDGET_H
