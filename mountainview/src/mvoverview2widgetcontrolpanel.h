#ifndef MVOVERVIEW2WIDGETCONTROLPANEL_H
#define MVOVERVIEW2WIDGETCONTROLPANEL_H

#include <QWidget>


class MVOverview2WidgetControlPanelPrivate;
class MVOverview2WidgetControlPanel : public QWidget
{
	Q_OBJECT
public:
	friend class MVOverview2WidgetControlPanelPrivate;
	MVOverview2WidgetControlPanel(QWidget *parent=0);
	virtual ~MVOverview2WidgetControlPanel();

	QVariant getParameterValue(QString name);
signals:
	void signalButtonClicked(QString str);
private slots:
	void slot_button_clicked();
private:
	MVOverview2WidgetControlPanelPrivate *d;
};

#endif // MVOVERVIEW2WIDGETCONTROLPANEL_H
