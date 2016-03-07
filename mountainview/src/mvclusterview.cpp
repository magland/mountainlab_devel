#include "mvclusterview.h"
#include <QImage>
#include <QColor>
#include <QPainter>
#include <stdio.h>
#include <QMouseEvent>
#include "affinetransformation.h"
#include <math.h>

class MVClusterViewPrivate {
public:
	MVClusterView *q;
	Mda m_data;
	Mda m_data_proj;
	bool m_data_proj_needed;
	Mda m_data_trans;
	bool m_data_trans_needed;

	QImage m_grid_image;
	Mda m_heat_map_grid;
	Mda m_point_grid;
	int m_grid_N1,m_grid_N2;
	double m_grid_delta;
	bool m_grid_update_needed;
	QPointF m_anchor_point;
	AffineTransformation m_anchor_transformation;

	Mda proj_matrix; //3xnum_features
	AffineTransformation m_transformation; //3x4

	void compute_data_proj();
	void compute_data_trans();
	void update_grid();
	void coord2gridindex(double x0,double y0,double &i1,double &i2);
	QColor get_heat_map_color(double val);
};

MVClusterView::MVClusterView(QWidget *parent) : QWidget(parent)
{
	d=new MVClusterViewPrivate;
	d->m_grid_N1=d->m_grid_N2=300;
	d->m_grid_delta=0.2;
	d->m_grid_update_needed=false;
	d->m_data_proj_needed=true;
	d->m_data_trans_needed=true;
	d->m_anchor_point=QPointF(-1,-1);
	d->m_transformation.setIdentity();
	this->setMouseTracking(true);
}

MVClusterView::~MVClusterView()
{
	delete d;
}

void MVClusterView::setData(const Mda &X)
{
	d->m_data=X;
	d->m_data_proj_needed=true;
	d->m_data_trans_needed=true;
	d->m_grid_update_needed=true;
	update();
}

QRectF compute_centered_square(QRectF R) {
	int margin=15;
	int W0=R.width()-margin*2;
	int H0=R.height()-margin*2;
	if (W0>H0) {
		return QRectF(margin+(W0-H0)/2,margin,H0,H0);
	}
	else {
		return QRectF(margin,margin+(H0-W0)/2,W0,W0);
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
	painter.fillRect(0,0,width(),height(),QColor(40,40,40));
	QRectF target=compute_centered_square(QRectF(0,0,width(),height()));
	painter.drawImage(target,d->m_grid_image);
	//QPen pen; pen.setColor(Qt::yellow);
	//painter.setPen(pen);
	//painter.drawRect(target);
}

void MVClusterView::mouseMoveEvent(QMouseEvent *evt)
{
	QPointF pt=evt->pos();
	if (d->m_anchor_point.x()>=0) {
		//We have the mouse button down!
		QPointF diff=pt-d->m_anchor_point;
		double factor=1.0/2;
		double deg_x=-diff.x()*factor;
		double deg_y=diff.y()*factor;
		d->m_transformation=d->m_anchor_transformation;
		d->m_transformation.rotateY(deg_x*M_PI/180);
		d->m_transformation.rotateX(deg_y*M_PI/180);
		d->m_data_trans_needed=true;
		d->m_grid_update_needed=true;
		update();
	}

}

void MVClusterView::mousePressEvent(QMouseEvent *evt)
{
	QPointF pt=evt->pos();
	d->m_anchor_point=pt;
	d->m_anchor_transformation=d->m_transformation;
}

void MVClusterView::mouseReleaseEvent(QMouseEvent *evt)
{
	Q_UNUSED(evt)
	d->m_anchor_point=QPointF(-1,-1);
}

void MVClusterView::wheelEvent(QWheelEvent *evt)
{
	double delta=evt->delta();
	double factor=1;
	if (delta>0) {
		factor=1.1;
	}
	else if (delta<0) {
		factor=1/1.1;
	}
	if (delta!=1) {
		d->m_transformation.scale(factor,factor,factor);
		d->m_data_trans_needed=true;
		d->m_grid_update_needed=true;
		update();
	}
}


void MVClusterViewPrivate::compute_data_proj()
{
	m_data_proj.allocate(3,m_data.N2());
	for (int i=0; i<m_data.N2(); i++) {
		m_data_proj.setValue(m_data.value(0,i),0,i);
		m_data_proj.setValue(m_data.value(1,i),1,i);
		m_data_proj.setValue(m_data.value(2,i),2,i);
	}
}

void MVClusterViewPrivate::compute_data_trans()
{
	int N2=m_data_proj.N2();
	m_data_trans.allocate(3,N2);
	double *AA=m_data_proj.dataPtr();
	double *BB=m_data_trans.dataPtr();
	double MM[16];
	m_transformation.getMatrixData(MM);
	int aaa=0;
	for (int i=0; i<N2; i++) {
		BB[aaa+0]=AA[aaa+0]*MM[0]+AA[aaa+1]*MM[1]+AA[aaa+2]*MM[2]  +  MM[3];
		BB[aaa+1]=AA[aaa+0]*MM[4]+AA[aaa+1]*MM[5]+AA[aaa+2]*MM[6]  +  MM[7];
		BB[aaa+2]=AA[aaa+0]*MM[8]+AA[aaa+1]*MM[9]+AA[aaa+2]*MM[10] +  MM[11];
		aaa+=3;
	}
	m_data_proj.write("/tmp/tmp_proj.mda");
	m_data_trans.write("/tmp/tmp_trans.mda");
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
	if (m_data_trans_needed) {
		compute_data_trans();
		m_data_trans_needed=false;
	}
	int kernel_rad=10;
	double kernel_tau=3;
	double kernel[(kernel_rad*2+1)*(kernel_rad*2+1)];
	{
		int aa=0;
		for (int dy=-kernel_rad; dy<=kernel_rad; dy++) {
			for (int dx=-kernel_rad; dx<=kernel_rad; dx++) {
				kernel[aa]=exp(-0.5*(dx*dx+dy*dy)/(kernel_tau*kernel_tau));
				aa++;
			}
		}
	}
	int N1=m_grid_N1,N2=m_grid_N2;
	m_point_grid.allocate(N1,N2);
	m_heat_map_grid.allocate(N1,N2);
	double *m_point_grid_ptr=m_point_grid.dataPtr();
	double *m_heat_map_grid_ptr=m_heat_map_grid.dataPtr();
	for (int j=0; j<m_data_trans.N2(); j++) {
		double x0=m_data_trans.value(0,j);
		double y0=m_data_trans.value(1,j);
		double factor=1;
		//if (m_data_trans.value(2,j)>0) {
			//factor=0.2;
		//}
		double i1,i2;
		coord2gridindex(x0,y0,i1,i2);
		int ii1=(int)(i1+0.5);
		int ii2=(int)(i2+0.5);
		if ((ii1-kernel_rad>=0)&&(ii1+kernel_rad<N1)&&(ii2-kernel_rad>=0)&&(ii2+kernel_rad<N2)) {
			m_point_grid_ptr[ii1+N1*ii2]=1;
			int aa=0;
			for (int dy=-kernel_rad; dy<=kernel_rad; dy++) {
				int ii2b=ii2+dy;
				for (int dx=-kernel_rad; dx<=kernel_rad; dx++) {
					int ii1b=ii1+dx;
					m_heat_map_grid_ptr[ii1b+N1*ii2b]+=kernel[aa]*factor;
					aa++;
				}
			}
		}
	}
	double max_heat_map_grid_val=compute_max(N1*N2,m_heat_map_grid.dataPtr());
	m_grid_image=QImage(N1,N2,QImage::Format_ARGB32);
	QColor white(255,255,255);
	QColor dark(50,50,50);
	for (int i2=0; i2<N2; i2++) {
		for (int i1=0; i1<N1; i1++) {
			double val=m_point_grid.value(i1,i2);
			if (val) {
				double val2=m_heat_map_grid.value(i1,i2)/max_heat_map_grid_val;
				QColor CC=get_heat_map_color(val2);
				//m_grid_image.setPixel(i1,i2,white.rgb());
				m_grid_image.setPixel(i1,i2,CC.rgb());
			}
			else {
				m_grid_image.setPixel(i1,i2,dark.rgb());
			}
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

QColor MVClusterViewPrivate::get_heat_map_color(double val)
{
	double r=0,g=0,b=0;
	if (val<0.2) {
		double tmp=(val-0)/0.2;
		r=200*(1-tmp)+150*tmp;
		b=200*(1-tmp)+255*tmp;
		g=0*(1-tmp)+0*tmp;
	}
	else if (val<0.4) {
		double tmp=(val-0.2)/0.2;
		r=150*(1-tmp)+0*tmp;
		b=255*(1-tmp)+255*tmp;
		g=0*(1-tmp)+100*tmp;
	}
	else if (val<0.6) {
		double tmp=(val-0.4)/0.2;
		r=0*(1-tmp)+255*tmp;
		b=255*(1-tmp)+0*tmp;
		g=100*(1-tmp)+20*tmp;
	}
	else if (val<0.8) {
		double tmp=(val-0.6)/0.2;
		r=255*(1-tmp)+255*tmp;
		b=0*(1-tmp)+0*tmp;
		g=20*(1-tmp)+255*tmp;
	}
	else if (val<=1.0) {
		double tmp=(val-0.8)/0.2;
		r=255*(1-tmp)+255*tmp;
		b=0*(1-tmp)+255*tmp;
		g=255*(1-tmp)+255*tmp;
	}

	return QColor((int)r,(int)g,(int)b);
}
