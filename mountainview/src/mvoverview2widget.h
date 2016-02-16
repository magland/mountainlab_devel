#ifndef MVOVERVIEW2WIDGET_H
#define MVOVERVIEW2WIDGET_H

#include <QWidget>
#include "mda.h"

class MVOverview2WidgetPrivate;
class MVOverview2Widget : public QWidget
{
	Q_OBJECT
public:
	friend class MVOverview2WidgetPrivate;
	MVOverview2Widget(QWidget *parent=0);
	virtual ~MVOverview2Widget();
	void setRawPath(const QString &path);
	void setFiringsPath(const QString &firings);
	void updateWidgets();

protected:
	void resizeEvent(QResizeEvent *evt);

signals:

public slots:

private slots:
	void slot_control_panel_button_clicked(QString str);

private:
	MVOverview2WidgetPrivate *d;
};

#endif // MVOVERVIEW2WIDGET_H
