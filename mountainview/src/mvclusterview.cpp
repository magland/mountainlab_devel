#include "mvclusterview.h"
#include <QImage>
#include <QColor>
#include <QPainter>

class MVClusterViewPrivate {
public:
	MVClusterView *q;
	Mda m_data;
	Mda m_data_proj;
	bool m_data_proj_needed;

	QImage m_grid_image;
	Mda m_grid;
	int m_grid_N1,m_grid_N2;
	double m_grid_delta;
	bool m_grid_update_needed;

	void compute_data_proj();
	void update_grid();
	void coord2gridindex(double x0,double y0,double &i1,double &i2);
};

MVClusterView::MVClusterView(QWidget *parent) : QWidget(parent)
{
	d=new MVClusterViewPrivate;
	d->m_grid_N1=d->m_grid_N2=300;
	d->m_grid_delta=0.1;
	d->m_grid_update_needed=false;
	d->m_data_proj_needed=true;
}

MVClusterView::~MVClusterView()
{
	delete d;
}

void MVClusterView::setData(const Mda &X)
{
	d->m_data=X;
	d->m_data_proj_needed=true;
	d->m_grid_update_needed=true;
	update();
}

QRectF compute_centered_square(QRectF R) {
	if (R.width()>R.height()) {
		return QRectF((R.width()-R.height())/2,0,R.height(),R.height());
	}
	else {
		return QRectF(0,(R.height()-R.width())/2,R.width(),R.width());
	}
}

void MVClusterView::paintEvent(QPaintEvent *evt)
{
	Q_UNUSED(evt)
	if (d->m_grid_update_needed) {
		d->update_grid();
		d->m_grid_update_needed=false;
	}

	QPainter painter(this);
	painter.fillRect(0,0,width(),height(),Qt::black);
	QRectF target=compute_centered_square(QRectF(0,0,width(),height()));
	painter.drawImage(target,d->m_grid_image);
	QPen pen; pen.setColor(Qt::yellow);
	painter.setPen(pen);
	painter.drawRect(target);
}


void MVClusterViewPrivate::compute_data_proj()
{
	m_data_proj.allocate(2,m_data.N2());
	for (int i=0; i<m_data.N2(); i++) {
		m_data_proj.setValue(m_data.value(0,i),0,i);
		m_data_proj.setValue(m_data.value(1,i),1,i);
	}
}

double compute_max(int N,double *X) {
	if (N==0) return 0;
	double ret=X[0];
	for (int i=0; i<N; i++) {
		if (X[i]>ret) ret=X[i];
	}
	return ret;
}

QColor make_color(double r,double g,double b) {
	if (r<0) r=0; if (r>1) r=1;
	if (g<0) g=0; if (g>1) g=1;
	if (b<0) b=0; if (b>1) b=1;
	return QColor((int)(r*255),(int)(g*255),(int)(b*255));
}

void MVClusterViewPrivate::update_grid()
{
	if (m_data_proj_needed) {
		compute_data_proj();
		m_data_proj_needed=false;
	}
	int N1=m_grid_N1,N2=m_grid_N2;
	m_grid.allocate(N1,N2);
	double *m_grid_ptr=m_grid.dataPtr();
	for (int j=0; j<m_data_proj.N2(); j++) {
		double x0=m_data_proj.value(0,j);
		double y0=m_data_proj.value(1,j);
		double i1,i2;
		coord2gridindex(x0,y0,i1,i2);
		int ii1=(int)(i1+0.5);
		int ii2=(int)(i2+0.5);
		if ((ii1>=0)&&(ii1<N1)&&(ii2>=0)&&(ii2<N2)) {
			m_grid_ptr[ii1+N1*ii2]++;
		}
	}
	double max_grid_val=compute_max(N1*N2,m_grid.dataPtr());
	m_grid_image=QImage(N1,N2,QImage::Format_ARGB32);
	for (int i2=0; i2<N2; i2++) {
		for (int i1=0; i1<N1; i1++) {
			double val=m_grid.value(i1,i2);
			val/=max_grid_val;
			double rr=val;
			double gg=val;
			double bb=val;
			QColor CC=make_color(rr,gg,bb);
			m_grid_image.setPixel(i1,i2,CC.rgb());
		}
	}
}

void MVClusterViewPrivate::coord2gridindex(double x0, double y0, double &i1, double &i2)
{
	int N1=m_grid_N1; int N2=m_grid_N2;
	int N1mid=(N1+1)/2-1; int N2mid=(N2+1)/2-1;
	i1=N1mid+x0/m_grid_delta;
	i2=N2mid+y0/m_grid_delta;
}
