#ifndef HISTOGRAMVIEW_H
#define HISTOGRAMVIEW_H

#include <QWidget>

class HistogramViewPrivate;
class HistogramView : public QWidget
{
	Q_OBJECT
public:
	friend class HistogramViewPrivate;
	explicit HistogramView(QWidget *parent = 0);
	virtual ~HistogramView();

    void setData(const QList<float> values);
	void setData(int N,float *values);
	void setBins(float bin_min,float bin_max,int num_bins);
	void autoSetBins(int num_bins);
	void setFillColor(const QColor &col);
	void setLineColor(const QColor &col);
    void setTitle(const QString &title);

protected:
	void paintEvent(QPaintEvent *evt);
    void mousePressEvent(QMouseEvent *evt);
	void mouseMoveEvent(QMouseEvent *evt);
    void enterEvent(QEvent *evt);
    void leaveEvent(QEvent *evt);
signals:
    void clicked();

public slots:

private:
	HistogramViewPrivate *d;
};

#endif // HISTOGRAMVIEW_H