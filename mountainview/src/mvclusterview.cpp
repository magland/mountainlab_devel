#include "mvclusterview.h"
#include <QImage>
#include <QColor>
#include <QPainter>
#include <stdio.h>
#include <QMouseEvent>
#include "affinetransformation.h"
#include <QTimer>
#include <math.h>

class MVClusterViewPrivate {
public:
	MVClusterView *q;
	Mda m_data;
	Mda m_data_proj;
	bool m_data_proj_needed;
	Mda m_data_trans;
	bool m_data_trans_needed;
	int m_current_event_index;
	int m_mode;

	QImage m_grid_image;
	QRectF m_image_target;
	Mda m_heat_map_grid;
	Mda m_point_grid;
	int m_grid_N1,m_grid_N2;
	double m_grid_delta;
	bool m_grid_update_needed;
	QPointF m_anchor_point;
	QPointF m_last_mouse_release_point;
	AffineTransformation m_anchor_transformation;
	bool m_moved_from_anchor;
	QList<double> m_times;
	QList<int> m_labels;
	QSet<int> m_closest_inds_to_exclude;
	QList<QColor> m_label_colors;
	bool m_emit_transformation_changed_scheduled;

	Mda proj_matrix; //3xnum_features
	AffineTransformation m_transformation; //3x4

	void compute_data_proj();
	void compute_data_trans();
	void update_grid();
	void coord2gridindex(double x0,double y0,double &i1,double &i2);
	QPointF pixel2coord(QPointF pix);
	QPointF coord2pixel(QPointF coord);
	QColor get_heat_map_color(double val);
	QColor get_label_color(int label);
	int find_closest_event_index(double x,double y,const QSet<int> &inds_to_exclude);
	void set_current_event_index(int ind,bool do_emit=true);
	void schedule_emit_transformation_changed();
};

MVClusterView::MVClusterView(QWidget *parent) : QWidget(parent)
{
	d=new MVClusterViewPrivate;
	d->q=this;
	d->m_grid_N1=d->m_grid_N2=300;
	d->m_grid_delta=8.0/d->m_grid_N1;
	d->m_grid_update_needed=false;
	d->m_data_proj_needed=true;
	d->m_data_trans_needed=true;
	d->m_anchor_point=QPointF(-1,-1);
	d->m_last_mouse_release_point=QPointF(-1,-1);
	d->m_transformation.setIdentity();
	d->m_moved_from_anchor=false;
	d->m_current_event_index=-1;
	d->m_mode=MVCV_MODE_LABEL_COLORS;
	d->m_emit_transformation_changed_scheduled=false;
	this->setMouseTracking(true);

	QList<QString> color_strings;
	color_strings << "#F7977A" << "#FDC68A"
				  << "#C4DF9B" << "#82CA9D"
				  << "#6ECFF6" << "#8493CA"
				  << "#A187BE" << "#F49AC2"
				  << "#F9AD81" << "#FFF79A"
				  << "#A2D39C" << "#7BCDC8"
				  << "#7EA7D8" << "#8882BE"
				  << "#BC8DBF" << "#F6989D";
	for (int i=0; i<color_strings.size(); i++) {
		d->m_label_colors << QColor(color_strings[i]);
	}
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

void MVClusterView::setTimes(const QList<double> &times)
{
	d->m_times=times;
}

void MVClusterView::setLabels(const QList<int> &labels)
{
	d->m_labels=labels;
}

void MVClusterView::setMode(int mode)
{
	d->m_mode=mode;
	update();
}

void MVClusterView::setCurrentEvent(MVEvent evt,bool do_emit)
{
	if ((d->m_times.isEmpty())||(d->m_labels.isEmpty())) return;
	for (int i=0; i<d->m_times.count(); i++) {
		if ((d->m_times[i]==evt.time)&&(d->m_labels.value(i)==evt.label)) {
			d->set_current_event_index(i,do_emit);
			return;
		}
	}
	d->set_current_event_index(-1);
}

MVEvent MVClusterView::currentEvent()
{
	MVEvent ret;
	if (d->m_current_event_index<0) {
		ret.time=-1;
		ret.label=-1;
		return ret;
	}
	ret.time=d->m_times.value(d->m_current_event_index);
	ret.label=d->m_labels.value(d->m_current_event_index);
	return ret;
}

int MVClusterView::currentEventIndex()
{
	return d->m_current_event_index;
}

AffineTransformation MVClusterView::transformation()
{
	return d->m_transformation;
}

void MVClusterView::setTransformation(const AffineTransformation &T)
{
	if (d->m_transformation.equals(T)) return;
	d->m_transformation=T;
	d->m_data_trans_needed=true;
	d->m_grid_update_needed=true;
	update();
	//do not emit to avoid excessive signals
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
	d->m_image_target=target;
	//QPen pen; pen.setColor(Qt::yellow);
	//painter.setPen(pen);
	//painter.drawRect(target);

	if (d->m_current_event_index>=0) {
		double x=d->m_data_trans.value(0,d->m_current_event_index);
		double y=d->m_data_trans.value(1,d->m_current_event_index);
		QPointF pix=d->coord2pixel(QPointF(x,y));
		painter.setBrush(QBrush(Qt::darkGreen));
		painter.drawEllipse(pix,6,6);
	}
}

void MVClusterView::mouseMoveEvent(QMouseEvent *evt)
{
	QPointF pt=evt->pos();
	if (d->m_anchor_point.x()>=0) {
		//We have the mouse button down!
		QPointF diff=pt-d->m_anchor_point;
		if ((qAbs(diff.x())>=5)||(qAbs(diff.y())>=5)||(d->m_moved_from_anchor)) {
			d->m_moved_from_anchor=true;
			double factor=1.0/2;
			double deg_x=-diff.x()*factor;
			double deg_y=diff.y()*factor;
			d->m_transformation=d->m_anchor_transformation;
			d->m_transformation.rotateY(deg_x*M_PI/180);
			d->m_transformation.rotateX(deg_y*M_PI/180);
			d->schedule_emit_transformation_changed();
			d->m_data_trans_needed=true;
			d->m_grid_update_needed=true;
			update();
		}
	}

}

void MVClusterView::mousePressEvent(QMouseEvent *evt)
{
	QPointF pt=evt->pos();
	d->m_anchor_point=pt;
	d->m_anchor_transformation=d->m_transformation;
	d->m_moved_from_anchor=false;
}

void MVClusterView::mouseReleaseEvent(QMouseEvent *evt)
{
	Q_UNUSED(evt)
	d->m_anchor_point=QPointF(-1,-1);
	if (evt->pos()!=d->m_last_mouse_release_point) d->m_closest_inds_to_exclude.clear();
	if (!d->m_moved_from_anchor) {
		QPointF coord=d->pixel2coord(evt->pos());
		int ind=d->find_closest_event_index(coord.x(),coord.y(),d->m_closest_inds_to_exclude);
		if (ind>=0) d->m_closest_inds_to_exclude.insert(ind);
		d->set_current_event_index(ind);
	}
	d->m_last_mouse_release_point=evt->pos();
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
		d->schedule_emit_transformation_changed();
		update();
	}
}

void MVClusterView::slot_emit_transformation_changed()
{
	d->m_emit_transformation_changed_scheduled=false;
	emit transformationChanged();
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
	double *m_point_grid_ptr=m_point_grid.dataPtr();

	if (m_mode==MVCV_MODE_HEAT_DENSITY) {
		m_heat_map_grid.allocate(N1,N2);
	}
	double *m_heat_map_grid_ptr=m_heat_map_grid.dataPtr();

	Mda z_grid;
	if (m_mode==MVCV_MODE_LABEL_COLORS) {
		z_grid.allocate(N1,N2);
	}
	double *z_grid_ptr=z_grid.dataPtr();

	for (int j=0; j<m_data_trans.N2(); j++) {
		double x0=m_data_trans.value(0,j);
		double y0=m_data_trans.value(1,j);
		double z0=m_data_trans.value(2,j);
		double factor=1;
		//if (m_data_trans.value(2,j)>0) {
			//factor=0.2;
		//}
		double i1,i2;
		coord2gridindex(x0,y0,i1,i2);
		int ii1=(int)(i1+0.5);
		int ii2=(int)(i2+0.5);
		if ((ii1-kernel_rad>=0)&&(ii1+kernel_rad<N1)&&(ii2-kernel_rad>=0)&&(ii2+kernel_rad<N2)) {
			int iiii=ii1+N1*ii2;
			if (m_mode==MVCV_MODE_LABEL_COLORS) {
				if ((m_point_grid_ptr[iiii]==0)||(z_grid_ptr[iiii]>z0)) {
					m_point_grid_ptr[iiii]=m_labels.value(j);
					z_grid_ptr[iiii]=z0;
				}
			}
			else {
				m_point_grid_ptr[iiii]=1;
			}

			if (m_mode==MVCV_MODE_HEAT_DENSITY) {
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
	}

	m_grid_image=QImage(N1,N2,QImage::Format_ARGB32);
	QColor dark(50,50,50);

	if (m_mode==MVCV_MODE_HEAT_DENSITY) {
		double max_heat_map_grid_val=compute_max(N1*N2,m_heat_map_grid.dataPtr());
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
	else {
		for (int i2=0; i2<N2; i2++) {
			for (int i1=0; i1<N1; i1++) {
				double val=m_point_grid.value(i1,i2);
				if (val) {
					QColor CC=get_label_color((int)val);
					m_grid_image.setPixel(i1,i2,CC.rgb());
				}
				else {
					m_grid_image.setPixel(i1,i2,dark.rgb());
				}
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

QPointF MVClusterViewPrivate::pixel2coord(QPointF pt)
{
	int N1=m_grid_N1; int N2=m_grid_N2;
	int N1mid=(N1+1)/2-1; int N2mid=(N2+1)/2-1;
	double pctx=(pt.x()-m_image_target.x())/(m_image_target.width());
	double pcty=(pt.y()-m_image_target.y())/(m_image_target.height());
	double xx=(-N1mid+pctx*N1)*m_grid_delta;
	double yy=(-N2mid+pcty*N2)*m_grid_delta;
	return QPointF(xx,yy);
}

QPointF MVClusterViewPrivate::coord2pixel(QPointF coord)
{
	int N1=m_grid_N1; int N2=m_grid_N2;
	int N1mid=(N1+1)/2-1; int N2mid=(N2+1)/2-1;
	double xx=coord.x();
	double yy=coord.y();
	double pctx=(xx/m_grid_delta+N1mid)/N1;
	double pcty=(yy/m_grid_delta+N2mid)/N2;
	double pt_x=pctx*m_image_target.width()+m_image_target.x();
	double pt_y=pcty*m_image_target.height()+m_image_target.y();
	return QPointF(pt_x,pt_y);
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

QColor MVClusterViewPrivate::get_label_color(int label)
{
	return m_label_colors[label % m_label_colors.size()];
}

int MVClusterViewPrivate::find_closest_event_index(double x, double y,const QSet<int> &inds_to_exclude)
{
	double best_distsqr=0;
	int best_ind=0;
	for (int i=0; i<m_data_trans.N2(); i++) {
		if (!inds_to_exclude.contains(i)) {
			double distx=m_data_trans.value(0,i)-x;
			double disty=m_data_trans.value(1,i)-y;
			double distsqr=distx*distx+disty*disty;
			if ((distsqr<best_distsqr)||(i==0)) {
				best_distsqr=distsqr;
				best_ind=i;
			}
		}
	}
	return best_ind;
}

void MVClusterViewPrivate::set_current_event_index(int ind,bool do_emit)
{
	if (m_current_event_index==ind) return;
	m_current_event_index=ind;
	if (do_emit) {
		emit q->currentEventChanged();
	}
	q->update();
}

void MVClusterViewPrivate::schedule_emit_transformation_changed()
{
	if (m_emit_transformation_changed_scheduled) return;
	m_emit_transformation_changed_scheduled=true;
	QTimer::singleShot(100,q,SLOT(slot_emit_transformation_changed()));
}

