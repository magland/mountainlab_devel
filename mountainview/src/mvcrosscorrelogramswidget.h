#ifndef MVCROSSCORRELOGRAMSWIDGET_H
#define MVCROSSCORRELOGRAMSWIDGET_H

#include <QWidget>

class MVCrossCorrelogramsWidgetPrivate;
class MVCrossCorrelogramsWidget : public QWidget
{
	Q_OBJECT
public:
	friend class MVCrossCorrelogramsWidgetPrivate;
	MVCrossCorrelogramsWidget();
	virtual ~MVCrossCorrelogramsWidget();

	void setCrossCorrelogramsPath(const QString &path);
	void updateWidget();

	int currentUnit();
	void setCurrentUnit(int num);

signals:
	void currentUnitChanged();

private slots:
	void slot_histogram_view_clicked();

private:
	MVCrossCorrelogramsWidgetPrivate *d;
};

#endif // MVCROSSCORRELOGRAMSWIDGET_H
