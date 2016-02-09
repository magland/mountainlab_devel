#ifndef MVMERGEWIDGET_H
#define MVMERGEWIDGET_H

#include <QWidget>
#include "mda.h"

class MVMergeWidgetPrivate;
class MVMergeWidget : public QWidget
{
public:
	friend class MVMergeWidgetPrivate;
	MVMergeWidget(QWidget *parent=0);
	virtual ~MVMergeWidget();
	void setRawPath(const QString &path);
	void setClustersPath(const QString &clusters);
	void setCorrelationMatrix(const Mda &CM);

protected:
	void resizeEvent(QResizeEvent *evt);

signals:

public slots:

private:
	MVMergeWidgetPrivate *d;
};

#endif // MVMERGEWIDGET_H
