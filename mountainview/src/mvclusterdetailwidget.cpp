#include "mvclusterdetailwidget.h"

#include <QPainter>
#include "mvutils.h"
#include "msutils.h"
#include <QProgressDialog>
#include <QTime>
#include <QMap>
#include <QDebug>
#include <QMouseEvent>

struct ClusterData {
	int k;
	int channel;
	Mda template0;
	QList<int> inds;
	QList<long> times;
	QList<double> peaks;
};

struct ChannelSpacingInfo {
	QList<double> channel_locations;
	double vert_scaling_factor;
};

class ClusterView {
public:
	friend class MVClusterDetailWidgetPrivate;
	friend class MVClusterDetailWidget;
	ClusterView(MVClusterDetailWidget *q0,MVClusterDetailWidgetPrivate *d0) {q=q0; d=d0; m_T=1; m_CD=0; m_highlighted=false; m_hovered=false;}
	void setClusterData(ClusterData *CD) {m_CD=CD;}
	void setChannelSpacingInfo(const ChannelSpacingInfo &csi) {m_csi=csi; m_T=0;}
	void setHighlighted(bool val) {m_highlighted=val;}
	void setHovered(bool val) {m_hovered=val;}
	void paint(QPainter *painter,QRectF rect);
	double spaceNeeded();
	ClusterData *clusterData() {return m_CD;}
	QRectF rect() {return m_rect;}

	MVClusterDetailWidget *q;
	MVClusterDetailWidgetPrivate *d;
private:
	ClusterData *m_CD;
	ChannelSpacingInfo m_csi;
	int m_T;
	QRectF m_rect;
	QRectF m_top_rect;
	QRectF m_template_rect;
	QRectF m_bottom_rect;
	bool m_highlighted;
	bool m_hovered;

	QPointF template_coord2pix(int m,double t,double val);
	QColor get_firing_rate_text_color(double rate);
};

class MVClusterDetailWidgetPrivate {
public:
	MVClusterDetailWidget *q;

	DiskReadMda m_raw;
	DiskReadMda m_firings;
	double m_sampling_freq;

	bool m_calculations_needed;
	int m_clip_size;
	QList<ClusterData> m_cluster_data;

	double m_vscale_factor;

	QProgressDialog *m_progress_dialog;
	QMap<QString,QColor> m_colors;
	QList<QColor> m_channel_colors;
	double m_total_time_sec;
	int m_current_k;
	int m_hovered_k;

	QList<ClusterView *> m_views;

	void do_calculations();
	void set_progress(QString title, QString text, float frac);
	void compute_total_time();
	void set_current_k(int k);
	void set_hovered_k(int k);
	int find_view_at(QPoint pos);
};

MVClusterDetailWidget::MVClusterDetailWidget(QWidget *parent) : QWidget(parent)
{
	d=new MVClusterDetailWidgetPrivate;
	d->q=this;
	d->m_calculations_needed=true;
	d->m_clip_size=100;
	d->m_progress_dialog=0;
	d->m_vscale_factor=1;
	d->m_total_time_sec=1;
	d->m_sampling_freq=0;
	d->m_current_k=-1;
	d->m_hovered_k=-1;

	d->m_colors["background"]=QColor(240,240,240);
	d->m_colors["frame1"]=QColor(245,245,245);
	d->m_colors["info_text"]=QColor(80,80,80);
	d->m_colors["view_background"]=QColor(245,245,245);
	d->m_colors["view_background_highlighted"]=QColor(250,220,200);
	d->m_colors["view_background_hovered"]=QColor(240,245,240);
	d->m_channel_colors << Qt::black;

	this->setFocusPolicy(Qt::StrongFocus);
	this->setMouseTracking(true);
}

MVClusterDetailWidget::~MVClusterDetailWidget()
{
	qDeleteAll(d->m_views);
	delete d;
}

void MVClusterDetailWidget::setRaw(DiskReadMda &X)
{
	d->m_raw=X;
	d->m_calculations_needed=true;
	d->compute_total_time();
	this->update();
}

void MVClusterDetailWidget::setFirings(DiskReadMda &X)
{
	d->m_firings=X;
	d->m_calculations_needed=true;
	this->update();
}

void MVClusterDetailWidget::setSamplingFrequency(double freq)
{
	d->m_sampling_freq=freq;
	d->compute_total_time();
	this->update();
}

void MVClusterDetailWidget::setChannelColors(const QList<QColor> &colors)
{
	d->m_channel_colors=colors;
	update();
}

int MVClusterDetailWidget::currentK()
{
	return d->m_current_k;
}

void MVClusterDetailWidget::setCurrentK(int k)
{
	d->set_current_k(k);
}

ChannelSpacingInfo compute_channel_spacing_info(QList<ClusterData> &cdata,double vscale_factor) {
	ChannelSpacingInfo info;
	info.vert_scaling_factor=1;
	if (cdata.count()==0) return info;
	int M=cdata[0].template0.N1();
	int T=cdata[0].template0.N2();
	double minval=0,maxval=0;
	for (int i=0; i<cdata.count(); i++) {
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				double val=cdata[i].template0.value(m,t);
				if (val<minval) minval=val;
				if (val>maxval) maxval=val;
			}
		}
	}
	double y0=0.5/M;
	for (int m=0; m<M; m++) {
		info.channel_locations << y0;
		y0+=1.0/M;
	}
	double maxabsval=qMax(maxval,-minval);
	info.vert_scaling_factor=0.5/M/maxabsval/vscale_factor;
	return info;
}

void MVClusterDetailWidget::paintEvent(QPaintEvent *evt)
{
	Q_UNUSED(evt)

	if (d->m_calculations_needed) {
		d->do_calculations();
		d->m_calculations_needed=false;
	}

	QPainter painter(this);
	painter.fillRect(0,0,width(),height(),d->m_colors["background"]);

	qDeleteAll(d->m_views);
	d->m_views.clear();
	for (int i=0; i<d->m_cluster_data.count(); i++) {
		ClusterData *CD=&d->m_cluster_data[i];
		ClusterView *V=new ClusterView(this,d);
		V->setHighlighted(CD->k==d->m_current_k);
		V->setHovered(CD->k==d->m_hovered_k);
		V->setClusterData(CD);
		d->m_views << V;
	}

	double total_space_needed=0;
	for (int i=0; i<d->m_views.count(); i++) {
		total_space_needed+=d->m_views[i]->spaceNeeded();
	}
	double space_ratio=width()/total_space_needed;

	ChannelSpacingInfo csi=compute_channel_spacing_info(d->m_cluster_data,d->m_vscale_factor);

	float x0=0;
	for (int i=0; i<d->m_views.count(); i++) {
		ClusterView *V=d->m_views[i];
		QRectF rect(x0,0,V->spaceNeeded()*space_ratio,height());
		V->setChannelSpacingInfo(csi);
		V->paint(&painter,rect);
		x0+=V->spaceNeeded()*space_ratio;
	}
}

void MVClusterDetailWidget::keyPressEvent(QKeyEvent *evt)
{

}

void MVClusterDetailWidget::mousePressEvent(QMouseEvent *evt)
{
	QPoint pt=evt->pos();
	int view_index=d->find_view_at(pt);
	if (view_index>=0) {
		int k=d->m_views[view_index]->clusterData()->k;
		if (d->m_current_k==k) d->set_current_k(-1);
		else d->set_current_k(k);
	}
	else {
		d->set_current_k(-1);
	}
}

void MVClusterDetailWidget::mouseReleaseEvent(QMouseEvent *evt)
{

}

void MVClusterDetailWidget::mouseMoveEvent(QMouseEvent *evt)
{
	QPoint pt=evt->pos();
	int view_index=d->find_view_at(pt);
	if (view_index>=0) {
		d->set_hovered_k(d->m_views[view_index]->clusterData()->k);
	}
	else {
		d->set_hovered_k(-1);
	}

}


void MVClusterDetailWidgetPrivate::do_calculations()
{

	int M=m_raw.N1();
	int N=m_raw.N2();
	int L=m_firings.N2();
	int T=m_clip_size;
	QList<long> times;
	QList<int> channels,labels;
	QList<double> peaks;

	Q_UNUSED(M)
	Q_UNUSED(N)

	for (int i=0; i<L; i++) {
		times << (long)m_firings.value(1,i)-1; //convert to 0-based indexing
		channels << (int)m_firings.value(0,i)-1; //convert to 0-based indexing
		labels << (int)m_firings.value(2,i);
		peaks << m_firings.value(3,i);
	}

	m_cluster_data.clear();
	//if we clear the cluster data, we need to make sure we delete the views! because there are pointers involved!
	//therefore this do_calculations() should only be called shortly before the views are redefined
	qDeleteAll(m_views);
	m_views.clear();

	int K=0;
	for (int i=0; i<L; i++) if (labels[i]>K) K=labels[i];

	for (int k=1; k<=K; k++) {
		set_progress("Computing Cluster Data","Computing Cluster Data",k*1.0/K);
		ClusterData CD;
		CD.k=k;
		CD.channel=0;
		for (int i=0; i<L; i++) {
			if (labels[i]==k) {
				CD.inds << i;
				CD.times << times[i];
				CD.channel=channels[i];
				CD.peaks << peaks[i];
			}
		}
		Mda clips_k=extract_clips(m_raw,CD.times,T);
		CD.template0=compute_mean_clip(clips_k);
		m_cluster_data << CD;
	}
}

void MVClusterDetailWidgetPrivate::set_progress(QString title, QString text, float frac)
{
	if (!m_progress_dialog) {
		m_progress_dialog=new QProgressDialog;
		m_progress_dialog->setCancelButton(0);
	}
	static QTime *timer=0;
	if (!timer) {
		timer=new QTime;
		timer->start();
		m_progress_dialog->show();
		m_progress_dialog->repaint();
	}
	if (timer->elapsed()>500) {
		timer->restart();
		if (!m_progress_dialog->isVisible()) {
			m_progress_dialog->show();
		}
		m_progress_dialog->setLabelText(text);
		m_progress_dialog->setWindowTitle(title);
		m_progress_dialog->setValue((int)(frac*100));
		m_progress_dialog->repaint();
	}
	if (frac>=1) {
		delete m_progress_dialog;
		m_progress_dialog=0;
	}
}

void MVClusterDetailWidgetPrivate::compute_total_time()
{
	m_total_time_sec=m_raw.N2()/m_sampling_freq;
}

void MVClusterDetailWidgetPrivate::set_current_k(int k)
{
	if (k==m_current_k) return;
	m_current_k=k;
	q->update();
	emit q->signalCurrentKChanged();
}

void MVClusterDetailWidgetPrivate::set_hovered_k(int k)
{
	if (k==m_hovered_k) return;
	m_hovered_k=k;
	q->update();
}

int MVClusterDetailWidgetPrivate::find_view_at(QPoint pos)
{
	for (int i=0; i<m_views.count(); i++) {
		if (m_views[i]->rect().contains(pos)) return i;
	}
	return -1;
}


void ClusterView::paint(QPainter *painter, QRectF rect)
{
	int mm=2;
	QRectF rect2(rect.x()+mm,rect.y()+mm,rect.width()-mm*2,rect.height()-mm*2);
	painter->setClipRect(rect);

	QColor background_color=d->m_colors["view_background"];
	if (m_highlighted) background_color=d->m_colors["view_background_highlighted"];
	else if (m_hovered) background_color=d->m_colors["view_background_hovered"];
	painter->fillRect(rect2,background_color);

	QPen pen_frame; pen_frame.setWidth(1); pen_frame.setColor(d->m_colors["frame1"]);
	painter->setPen(pen_frame);
	painter->drawRect(rect2);

	Mda template0=m_CD->template0;
	int M=template0.N1();
	int T=template0.N2();
	m_T=T;

	int top_height=20,bottom_height=40;
	m_rect=rect2;
	m_top_rect=QRectF(m_rect.x(),m_rect.y(),m_rect.width(),top_height);
	m_template_rect=QRectF(m_rect.x(),m_rect.y()+top_height,m_rect.width(),m_rect.height()-bottom_height-top_height);
	m_bottom_rect=QRectF(m_rect.x(),m_rect.y()+m_rect.height()-bottom_height,m_rect.width(),bottom_height);

	QPen pen; pen.setWidth(1);
	for (int m=0; m<M; m++) {
		QColor col=d->m_channel_colors.value(m%d->m_channel_colors.count());
		pen.setColor(col);
		painter->setPen(pen);
		QPainterPath path;
		for (int t=0; t<T; t++) {
			QPointF pt=template_coord2pix(m,t,template0.value(m,t));
			if (t==0) path.moveTo(pt);
			else path.lineTo(pt);
		}
		painter->drawPath(path);
	}

	QFont font=painter->font();
	QString txt;
	QRectF RR;

	txt=QString("%1").arg(m_CD->k);
	font.setPixelSize(16); pen.setColor(Qt::darkBlue);
	painter->setFont(font); painter->setPen(pen);
	painter->drawText(m_top_rect,Qt::AlignCenter|Qt::AlignBottom,txt);

	font.setPixelSize(11);
	int text_height=13;

	RR=QRectF(m_bottom_rect.x(),m_bottom_rect.y()+m_bottom_rect.height()-text_height,m_bottom_rect.width(),text_height);
	txt=QString("%1 spikes").arg(m_CD->inds.count());
	pen.setColor(d->m_colors["info_text"]);
	painter->setFont(font); painter->setPen(pen);
	painter->drawText(RR,Qt::AlignCenter|Qt::AlignBottom,txt);

	RR=QRectF(m_bottom_rect.x(),m_bottom_rect.y()+m_bottom_rect.height()-text_height*2,m_bottom_rect.width(),text_height);
	double rate=m_CD->inds.count()*1.0/d->m_total_time_sec;
	pen.setColor(get_firing_rate_text_color(rate));
	txt=QString("%1 sp/sec").arg(QString::number(rate,'g',2));
	painter->setFont(font); painter->setPen(pen);
	painter->drawText(RR,Qt::AlignCenter|Qt::AlignBottom,txt);
}

double ClusterView::spaceNeeded()
{
	return 1;
}

QPointF ClusterView::template_coord2pix(int m, double t, double val)
{
	double pcty=m_csi.channel_locations.value(m)+val*m_csi.vert_scaling_factor;
	double pctx=0;
	if (m_T) pctx=(t+0.5)/m_T;
	int margx=4;
	int margy=4;
	float x0=m_template_rect.x()+margx+pctx*(m_template_rect.width()-margx*2);
	float y0=m_template_rect.y()+margy+pcty*(m_template_rect.height()-margy*2);
	return QPointF(x0,y0);
}

QColor ClusterView::get_firing_rate_text_color(double rate)
{
	if (rate<=0.1) return QColor(220,220,220);
	if (rate<=1) return QColor(150,150,150);
	if (rate<=10) return QColor(0,50,0);
	return QColor(50,0,0);
}

