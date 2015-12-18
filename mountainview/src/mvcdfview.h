#ifndef MVCDFVIEW_H
#define MVCDFVIEW_H
#include <QVector>
#include <QWidget>
#include "mda.h"

class MVCdfViewPrivate;
class MVCdfView : public QWidget
{
	Q_OBJECT
public:
	friend class MVCdfViewPrivate;
	MVCdfView();
	virtual ~MVCdfView();
	void setTimesLabels(const QVector<int> &times,const QVector<int> &labels);
	void setTimesLabels(const Mda &times,const Mda &labels);
	void setCurrentLabel(int val);
	int currentLabel();
protected:
	void paintEvent(QPaintEvent *evt);
	void resizeEvent(QResizeEvent *evt);
	void mousePressEvent(QMouseEvent *evt);
signals:

private slots:

private:
	MVCdfViewPrivate *d;
};

#endif // MVCDFVIEW_H
