#ifndef FTOPTIONSWIDGET_H
#define FTOPTIONSWIDGET_H

#include <QWidget>


class FTOptionsWidgetPrivate;

class FTOptionsWidget : public QWidget
{
	Q_OBJECT
public:
	friend class FTOptionsWidgetPrivate;
	explicit FTOptionsWidget(QWidget *parent = 0);
	~FTOptionsWidget();

	bool showChannelNumbers();
	bool autoSelectChannels();
	bool normalizeIntensity();
	float brightness();

	void setShowChannelNumbers(bool val);

private:
	FTOptionsWidgetPrivate *d;

signals:
	void signalOptionsChanged();

public slots:
};

#endif // FTOPTIONSWIDGET_H
