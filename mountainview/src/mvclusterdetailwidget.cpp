#include "mvclusterdetailwidget.h"

#include <QPainter>
#include "mvutils.h"
#include "msutils.h"
#include <QProgressDialog>
#include <QTime>
#include <QMap>
#include <QDebug>
#include <QMouseEvent>
#include <QSet>

struct Subcluster {
    QList<int> inds;
};

struct ClusterData {
	int k;
	int channel;
    Mda templates;
	QList<int> inds;
    QList<Subcluster> subclusters;
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
	ClusterView(MVClusterDetailWidget *q0,MVClusterDetailWidgetPrivate *d0) {q=q0; d=d0; m_T=1; m_CD=0; m_highlighted=false; m_hovered=false; x_position_before_scaling=0;}
	void setClusterData(ClusterData *CD) {m_CD=CD;}
	void setChannelSpacingInfo(const ChannelSpacingInfo &csi) {m_csi=csi; m_T=0;}
	void setHighlighted(bool val) {m_highlighted=val;}
	void setHovered(bool val) {m_hovered=val;}
	void setSelected(bool val) {m_selected=val;}
	void paint(QPainter *painter,QRectF rect);
	double spaceNeeded();
	ClusterData *clusterData() {return m_CD;}
	QRectF rect() {return m_rect;}

	MVClusterDetailWidget *q;
	MVClusterDetailWidgetPrivate *d;
	double x_position_before_scaling;
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
	bool m_selected;

    QPointF template_coord2pix(int template_num,int m,double t,double val);
	QColor get_firing_rate_text_color(double rate);
};

class MVClusterDetailWidgetPrivate {
public:
	MVClusterDetailWidget *q;

	DiskReadMda m_raw;
	DiskReadMda m_firings;
    QList<int> m_group_numbers;
	double m_sampling_freq;

	bool m_calculations_needed;
	int m_clip_size;
	QList<ClusterData> m_cluster_data;

	double m_vscale_factor;
	double m_space_ratio;
	double m_scroll_x;

	QProgressDialog *m_progress_dialog;
	QMap<QString,QColor> m_colors;
	QList<QColor> m_channel_colors;
	double m_total_time_sec;
	int m_current_k;
	QSet<int> m_selected_ks;
	int m_hovered_k;
	double m_anchor_x; double m_anchor_scroll_x;

	QList<ClusterView *> m_views;

	void do_calculations();
	void set_progress(QString title, QString text, float frac);
	void compute_total_time();
	void set_current_k(int k);
	void set_hovered_k(int k);
	int find_view_index_at(QPoint pos);
	ClusterView *find_view_for_k(int k);
	int find_view_index_for_k(int k);
	void ensure_view_visible(ClusterView *V);
	void zoom(double factor);
};

MVClusterDetailWidget::MVClusterDetailWidget(QWidget *parent) : QWidget(parent)
{
	d=new MVClusterDetailWidgetPrivate;
	d->q=this;
	d->m_calculations_needed=true;
	d->m_clip_size=100;
	d->m_progress_dialog=0;
	d->m_vscale_factor=2;
	d->m_total_time_sec=1;
	d->m_sampling_freq=0;
	d->m_current_k=-1;
	d->m_hovered_k=-1;
	d->m_space_ratio=50;
	d->m_scroll_x=0;
	d->m_anchor_x=-1;
	d->m_anchor_scroll_x=-1;

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

int compute_max_label(DiskReadMda &firings) {
    int K=0;
    for (int i=0; i<firings.N2(); i++) {
        int val=(int)firings.value(2,i);
        if (val>K) K=val;
    }
    return K;
}

void MVClusterDetailWidget::setFirings(DiskReadMda &X,const QList<int> &group_numbers_in)
{
    QList<int> group_numbers=group_numbers_in;
    if (group_numbers.isEmpty()) {
        int K=compute_max_label(X);
        group_numbers << 0;
        for (int k=1; k<=K; k++) group_numbers << k;
    }
	d->m_firings=X;
    d->m_group_numbers=group_numbers;
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

void MVClusterDetailWidget::setColors(const QMap<QString, QColor> &colors)
{
	d->m_colors=colors;
}

int MVClusterDetailWidget::currentK()
{
	return d->m_current_k;
}

QList<int> MVClusterDetailWidget::selectedKs()
{
	return QList<int>::fromSet(d->m_selected_ks);
}

void MVClusterDetailWidget::setCurrentK(int k)
{
	d->set_current_k(k);
}

void MVClusterDetailWidget::setSelectedKs(const QList<int> &ks)
{
	d->m_selected_ks=ks.toSet();
	update();
}

ChannelSpacingInfo compute_channel_spacing_info(QList<ClusterData> &cdata,double vscale_factor) {
	ChannelSpacingInfo info;
	info.vert_scaling_factor=1;
	if (cdata.count()==0) return info;
    int M=cdata[0].templates.N1();
    int T=cdata[0].templates.N2();
	double minval=0,maxval=0;
	for (int i=0; i<cdata.count(); i++) {
        int NT=cdata[i].templates.N3();
        for (int ii=0; ii<NT; ii++) {
            for (int t=0; t<T; t++) {
                for (int m=0; m<M; m++) {
                    double val=cdata[i].templates.value(m,t,ii);
                    if (val<minval) minval=val;
                    if (val>maxval) maxval=val;
                }
            }
        }
	}
	double y0=0.5/M;
	for (int m=0; m<M; m++) {
		info.channel_locations << y0;
		y0+=1.0/M;
	}
	double maxabsval=qMax(maxval,-minval);
	info.vert_scaling_factor=0.5/M/maxabsval*vscale_factor;
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
		V->setSelected(d->m_selected_ks.contains(CD->k));
		V->setHovered(CD->k==d->m_hovered_k);
		V->setClusterData(CD);
		d->m_views << V;
	}

	double total_space_needed=0;
	for (int i=0; i<d->m_views.count(); i++) {
		total_space_needed+=d->m_views[i]->spaceNeeded();
	}
	if (d->m_scroll_x<0) d->m_scroll_x=0;
	if (total_space_needed*d->m_space_ratio-d->m_scroll_x<this->width()) {
		d->m_scroll_x=total_space_needed*d->m_space_ratio-this->width();
		if (d->m_scroll_x<0) d->m_scroll_x=0;
	}
	if ((d->m_scroll_x==0)&&(total_space_needed*d->m_space_ratio<this->width())) {
		d->m_space_ratio=this->width()/total_space_needed;
		if (d->m_space_ratio>300) d->m_space_ratio=300;
	}


	ChannelSpacingInfo csi=compute_channel_spacing_info(d->m_cluster_data,d->m_vscale_factor);

	float x0_before_scaling=0;
	for (int i=0; i<d->m_views.count(); i++) {
		ClusterView *V=d->m_views[i];
		QRectF rect(x0_before_scaling*d->m_space_ratio-d->m_scroll_x,0,V->spaceNeeded()*d->m_space_ratio,height());
		V->setChannelSpacingInfo(csi);
		V->paint(&painter,rect);
		V->x_position_before_scaling=x0_before_scaling;
		x0_before_scaling+=V->spaceNeeded();
	}
}

void MVClusterDetailWidget::keyPressEvent(QKeyEvent *evt)
{
	double factor=1.15;
	if (evt->key()==Qt::Key_Up) {
		d->m_vscale_factor*=factor;
		update();
	}
	else if (evt->key()==Qt::Key_Down) {
		d->m_vscale_factor/=factor;
		update();
	}
	else if ((evt->key()==Qt::Key_Plus)||(evt->key()==Qt::Key_Equal)) {
		d->zoom(1.1);
	}
	else if (evt->key()==Qt::Key_Minus) {
		d->zoom(1/1.1);
	}
}

void MVClusterDetailWidget::mousePressEvent(QMouseEvent *evt)
{
	QPoint pt=evt->pos();
	d->m_anchor_x=pt.x();
	d->m_anchor_scroll_x=d->m_scroll_x;
}

void MVClusterDetailWidget::mouseReleaseEvent(QMouseEvent *evt)
{
	QPoint pt=evt->pos();

	if ((d->m_anchor_x>=0)&&(qAbs(pt.x()-d->m_anchor_x)>5)) {
		d->m_scroll_x=d->m_anchor_scroll_x-(pt.x()-d->m_anchor_x);
		d->m_anchor_x=-1;
		update();
		return;
	}
	d->m_anchor_x=-1;

	if (evt->modifiers()&Qt::ControlModifier) {
		int view_index=d->find_view_index_at(pt);
		if (view_index>=0) {
			int k=d->m_views[view_index]->clusterData()->k;
            if (d->m_current_k==k) {
                d->set_current_k(-1);
            }
			if (d->m_selected_ks.contains(k)) {
				d->m_selected_ks.remove(k);
				emit signalSelectedKsChanged();
				update();
			}
			else {
				d->m_selected_ks.insert(k);
				emit signalSelectedKsChanged();
				update();
			}
		}
	}
	else {
		int view_index=d->find_view_index_at(pt);
		if (view_index>=0) {
			int k=d->m_views[view_index]->clusterData()->k;
			if (d->m_current_k==k) {

			}
			else {
				d->set_current_k(k);
				d->m_selected_ks.clear();
				d->m_selected_ks.insert(k);
				emit signalSelectedKsChanged();
				update();
			}
		}
		else {
			d->set_current_k(-1);
			d->m_selected_ks.clear();
			emit signalSelectedKsChanged();
			update();
		}
	}
}

void MVClusterDetailWidget::mouseMoveEvent(QMouseEvent *evt)
{
	QPoint pt=evt->pos();
	if ((d->m_anchor_x>=0)&&(qAbs(pt.x()-d->m_anchor_x)>5)) {
		d->m_scroll_x=d->m_anchor_scroll_x-(pt.x()-d->m_anchor_x);
		update();
		return;
	}

	int view_index=d->find_view_index_at(pt);
	if (view_index>=0) {
		d->set_hovered_k(d->m_views[view_index]->clusterData()->k);
	}
	else {
		d->set_hovered_k(-1);
	}

}

void MVClusterDetailWidget::wheelEvent(QWheelEvent *evt)
{
	int delta=evt->delta();
	double factor=1;
	if (delta>0) factor=1.1;
	else factor=1/1.1;
	d->zoom(factor);
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
        CD.templates=compute_mean_clip(clips_k);
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
	ClusterView *V=find_view_for_k(k);
	if (V) ensure_view_visible(V);
	q->update();
	emit q->signalCurrentKChanged();
}

void MVClusterDetailWidgetPrivate::set_hovered_k(int k)
{
	if (k==m_hovered_k) return;
	m_hovered_k=k;
	q->update();
}

int MVClusterDetailWidgetPrivate::find_view_index_at(QPoint pos)
{
	for (int i=0; i<m_views.count(); i++) {
		if (m_views[i]->rect().contains(pos)) return i;
	}
	return -1;
}

ClusterView *MVClusterDetailWidgetPrivate::find_view_for_k(int k)
{
	for (int i=0; i<m_views.count(); i++) {
		if (m_views[i]->clusterData()->k==k) return m_views[i];
	}
	return 0;
}

int MVClusterDetailWidgetPrivate::find_view_index_for_k(int k)
{
	for (int i=0; i<m_views.count(); i++) {
		if (m_views[i]->clusterData()->k==k) return i;
	}
	return -1;
}

void MVClusterDetailWidgetPrivate::ensure_view_visible(ClusterView *V)
{
	double x0=V->x_position_before_scaling*m_space_ratio;
	if (x0<m_scroll_x) {
		m_scroll_x=x0-100;
		if (m_scroll_x<0) m_scroll_x=0;
	}
	else if (x0>m_scroll_x+q->width()) {
		m_scroll_x=x0-q->width()+100;
	}
}

void MVClusterDetailWidgetPrivate::zoom(double factor)
{
	if ((m_current_k>=0)&&(find_view_for_k(m_current_k))) {
		ClusterView *view=find_view_for_k(m_current_k);
		double current_screen_x=view->x_position_before_scaling*m_space_ratio-m_scroll_x;
		m_space_ratio*=factor;
		m_scroll_x=view->x_position_before_scaling*m_space_ratio-current_screen_x;
		if (m_scroll_x<0) m_scroll_x=0;
	}
	else {
		m_space_ratio*=factor;
	}
	q->update();
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

	QPen pen_frame; pen_frame.setWidth(1);
	pen_frame.setColor(d->m_colors["frame1"]);
	if (m_selected) pen_frame.setColor(d->m_colors["view_frame_selected"]);
	painter->setPen(pen_frame);
	painter->drawRect(rect2);

    Mda templates=m_CD->templates;
    int M=templates.N1();
    int T=templates.N2();
    int NT=templates.N3();
	m_T=T;

	int top_height=20,bottom_height=40;
	m_rect=rect;
	m_top_rect=QRectF(rect2.x(),rect2.y(),rect2.width(),top_height);
	m_template_rect=QRectF(rect2.x(),rect2.y()+top_height,rect2.width(),rect2.height()-bottom_height-top_height);
	m_bottom_rect=QRectF(rect2.x(),rect2.y()+rect2.height()-bottom_height,rect2.width(),bottom_height);

	QPen pen; pen.setWidth(1);
    for (int ii=0; ii<NT; ii++) {
        for (int m=0; m<M; m++) {
            QColor col=d->m_channel_colors.value(m%d->m_channel_colors.count());
            pen.setColor(col);
            painter->setPen(pen);
            QPainterPath path;
            for (int t=0; t<T; t++) {
                QPointF pt=template_coord2pix(ii,m,t,templates.value(m,t,ii));
                if (t==0) path.moveTo(pt);
                else path.lineTo(pt);
            }
            painter->drawPath(path);
        }
    }

	QFont font=painter->font();
	QString txt;
	QRectF RR;

	bool compressed_info=false;
	if (rect2.width()<60) compressed_info=true;

	txt=QString("%1").arg(m_CD->k);
	font.setPixelSize(16);
	if (compressed_info) font.setPixelSize(12);
	pen.setColor(Qt::darkBlue);
	painter->setFont(font); painter->setPen(pen);
	painter->drawText(m_top_rect,Qt::AlignCenter|Qt::AlignBottom,txt);

	font.setPixelSize(11);
	int text_height=13;

	if (!compressed_info) {
		RR=QRectF(m_bottom_rect.x(),m_bottom_rect.y()+m_bottom_rect.height()-text_height,m_bottom_rect.width(),text_height);
		txt=QString("%1 spikes").arg(m_CD->inds.count());
		pen.setColor(d->m_colors["info_text"]);
		painter->setFont(font); painter->setPen(pen);
		painter->drawText(RR,Qt::AlignCenter|Qt::AlignBottom,txt);
	}

	RR=QRectF(m_bottom_rect.x(),m_bottom_rect.y()+m_bottom_rect.height()-text_height*2,m_bottom_rect.width(),text_height);
	double rate=m_CD->inds.count()*1.0/d->m_total_time_sec;
	pen.setColor(get_firing_rate_text_color(rate));
	if (!compressed_info) txt=QString("%1 sp/sec").arg(QString::number(rate,'g',2));
	else txt=QString("%1").arg(QString::number(rate,'g',2));
	painter->setFont(font); painter->setPen(pen);
	painter->drawText(RR,Qt::AlignCenter|Qt::AlignBottom,txt);
}

double ClusterView::spaceNeeded()
{
	return 1;
}

QPointF ClusterView::template_coord2pix(int template_num,int m, double t, double val)
{
    int num_templates=this->m_CD->templates.N3();
    double pcty=m_csi.channel_locations.value(m)-val*m_csi.vert_scaling_factor; //negative because (0,0) is top-left, not bottom-right
	double pctx=0;
	if (m_T) pctx=(t+0.5)/m_T;
    pctx=template_num/num_templates+pctx/num_templates;
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


