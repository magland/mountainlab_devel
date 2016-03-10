#include "mvfiringrateview.h"

#include <QPainter>
#include <QDebug>
#include <math.h>
#include <QHBoxLayout>
#include <QTimer>
#include "sstimeseriesview.h"

class MVFiringRateViewPrivate {
public:
    MVFiringRateView *q;

    //set by user
    DiskReadMda m_firings;
    QList<int> m_visible_labels;
    double m_window_width_sec; //seconds
    double m_smoothing_kernel_size;
    double m_bin_size; //timepoints
    double m_sampling_freq;
    bool m_log_mode;

    //internal
    double m_max_timepoint;
    int m_K;
    Mda m_event_counts;
    double m_max_firing_rate;
    SSTimeSeriesView *m_view;
    bool m_update_scheduled;

    //settings
    double m_hmargin,m_vmargin;

    void do_compute();
    QPointF coord2pix(const QPointF &p);
    QColor get_label_color(int label);
    void schedule_update();
};

MVFiringRateView::MVFiringRateView()
{
    d=new MVFiringRateViewPrivate;
    d->q=this;
    d->m_window_width_sec=1; //seconds
    d->m_smoothing_kernel_size=5; //seconds
    d->m_sampling_freq=30000;
    d->m_bin_size=d->m_window_width_sec*d->m_sampling_freq; //timepoints
    d->m_max_timepoint=0;
    d->m_max_firing_rate=0;
    d->m_log_mode=true;
    d->m_update_scheduled=false;

    d->m_hmargin=25;
    d->m_vmargin=25;

    d->m_view=new SSTimeSeriesView;
    d->m_view->initialize();
    d->m_view->plot()->setFixedVerticalChannelSpacing(0);

    QHBoxLayout *hlayout=new QHBoxLayout;
    hlayout->addWidget(d->m_view);
    this->setLayout(hlayout);
}

MVFiringRateView::~MVFiringRateView()
{
    delete d;
}

void MVFiringRateView::setFirings(const DiskReadMda &X)
{
    d->m_firings=X;
    d->m_max_timepoint=0;
    d->m_K=0;
    for (int i=0; i<X.N2(); i++) {
        if (d->m_firings.value(1,i)>d->m_max_timepoint) d->m_max_timepoint=d->m_firings.value(1,i);
        if (d->m_firings.value(2,i)>d->m_K) d->m_K=(int)d->m_firings.value(2,i);
    }
    d->schedule_update();
}

void MVFiringRateView::setVisibleLabels(const QList<int> &X)
{
    d->m_visible_labels=X;
    d->schedule_update();
}

void MVFiringRateView::setWindowWidth(double sec)
{
    d->m_window_width_sec=sec;
    d->m_bin_size=d->m_window_width_sec*d->m_sampling_freq;
    d->schedule_update();
}

void MVFiringRateView::setSamplingFreq(double ff)
{
    d->m_sampling_freq=ff;
    d->m_bin_size=d->m_window_width_sec*d->m_sampling_freq;
    d->schedule_update();
}

void MVFiringRateView::setCurrentEvent(MVEvent evt)
{
    Q_UNUSED(evt)
    //finish this!
}

QList<double> do_smoothing(const QList<double> &X,double sigma) {
    int krad=(int)(sigma*4);
    double kernel[krad*2+1];
    //double kernel_sum=0;
    for (int dt=-krad; dt<=krad; dt++) {
        kernel[dt+krad]=exp(-0.5*dt*dt/(sigma*sigma));
        //kernel_sum+=kernel[dt+krad];
    }
    //for (int dt=-krad; dt<=krad; dt++) {
    //    kernel[dt+krad]/=kernel_sum;
    //}

    QList<double> Y;
    for (int t=0; t<X.count(); t++) {
        double val=0;
        double denom=0;
        for (int dt=-krad; dt<=krad; dt++) {
            if ((0<=t+dt)&&(t+dt<X.count())) {
                val+=X[t+dt]*kernel[dt+krad];
                denom+=kernel[dt+krad];
            }
        }
        Y << val/denom;
    }

    return Y;
}

void MVFiringRateView::slot_update()
{
    d->m_update_scheduled=false;
    d->do_compute();
    int Kv=d->m_event_counts.N1();
    long num_bins=d->m_event_counts.N2();
    Mda data; data.allocate(Kv,num_bins);
    for (int k=0; k<Kv; k++) {
        QList<double> rates;
        for (int i=0; i<num_bins; i++) {
            rates << d->m_event_counts.value(k,i)/d->m_window_width_sec;
        }
        QList<double> rates_smoothed=do_smoothing(rates,d->m_smoothing_kernel_size/d->m_window_width_sec);
        for (int i=0; i<rates_smoothed.count(); i++) {
            double val=0;
            if (rates_smoothed[i]) val=log(rates_smoothed[i]);
            data.setValue(val,k,i);
        }
    }
    DiskArrayModel *DAM=new DiskArrayModel();
    DAM->setFromMda(data);
    d->m_view->setData(DAM,this);
}



/*
void MVFiringRateView::paintEvent(QPaintEvent *evt)
{
    Q_UNUSED(evt);
    if (d->m_compute_needed) {
        d->do_compute();
        d->m_compute_needed=false;
    }

    QPainter painter(this);
    painter.fillRect(0,0,width(),height(),QColor(50,50,50));

    int Kv=d->m_event_counts.N1();
    long num_bins=d->m_event_counts.N2();
    for (int k=0; k<Kv; k++) {
        QPainterPath path;
        QList<double> rates;
        for (int i=0; i<num_bins; i++) {
            rates << d->m_event_counts.value(k,i)/d->m_window_width_sec;
        }
        QList<double> rates_smoothed=do_smoothing(rates,d->m_smoothing_kernel_size/d->m_window_width_sec);
        for (int i=0; i<rates_smoothed.count(); i++) {
            double rate=rates_smoothed[i];
            double time0=(i+0.5)*d->m_bin_size;
            QPointF pt=d->coord2pix(QPointF(time0,rate));
            if (i==0) path.moveTo(pt);
            else path.lineTo(pt);
        }
        QColor col=d->get_label_color(k);
        QPen pen; pen.setColor(col);
        painter.setPen(pen);
        painter.drawPath(path);
    }
}
*/


void MVFiringRateViewPrivate::do_compute()
{
    long num_bins=2+ (long)(m_max_timepoint/m_bin_size);
    int Kv=m_visible_labels.count();
    QList<int> label_map;
    for (int k=0; k<=m_K; k++) label_map << -1;
    for (int ii=0; ii<Kv; ii++) {
        if ((m_visible_labels[ii]>=0)&&(m_visible_labels[ii]<label_map.count())) {
            label_map[m_visible_labels[ii]]=ii;
        }
        else {
            qWarning() << "Unexpected problem" << __FILE__ << __FUNCTION__ << __LINE__;
        }

    }
    m_event_counts.allocate(Kv,num_bins);
    for (int i=0; i<m_firings.N2(); i++) {
        int k=(int)m_firings.value(2,i);
        int ii=label_map[k];
        if (ii>=0) {
            double t=m_firings.value(1,i);
            long bin=t/m_bin_size;
            m_event_counts.setValue(m_event_counts.value(ii,bin)+1,ii,bin);
        }
    }
    m_max_firing_rate=0;
    for (int i=0; i<num_bins; i++) {
        for (int j=0; j<Kv; j++) {
            double rate=m_event_counts.value(j,i)/m_window_width_sec;
            if (rate>m_max_firing_rate) m_max_firing_rate=rate;
        }
    }
}

QPointF MVFiringRateViewPrivate::coord2pix(const QPointF &p)
{
    if ((m_max_timepoint==0)||(m_max_firing_rate==0)) return QPointF(0,0);
    double pctx=p.x()/m_max_timepoint;
    double pcty;

    if (m_log_mode) {
        double v1=log(m_max_firing_rate/10000);
        double v2=log(m_max_firing_rate);
        if (p.y()) {
            pcty=(log(p.y())-v1)/(v2-v1);
        }
        else {
            pcty=0;
        }
    }
    else {
        pcty=p.y()/m_max_firing_rate;
    }
    if (pcty<0) pcty=0;
    if (pcty>1) pcty=1;

    double x0=m_hmargin+pctx*(q->width()-2*m_hmargin);
    double y0=m_vmargin+(1-pcty)*(q->height()-2*m_vmargin);

    return QPointF(x0,y0);
}

void MVFiringRateViewPrivate::schedule_update()
{
    if (m_update_scheduled) return;
    m_update_scheduled=true;
    QTimer::singleShot(500,q,SLOT(slot_update()));
}
