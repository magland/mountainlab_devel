#include "mvfiringrateview.h"

#include <QPainter>
#include <QDebug>
#include <math.h>
#include <QHBoxLayout>
#include <QTimer>
#include "sstimeseriesview.h"
#include "mvutils.h"
#include "msutils.h"

class MVFiringRateViewPrivate {
public:
    MVFiringRateView *q;

    //set by user
    Mda m_firings;
    double m_sampling_freq;
    QList<Epoch> m_epochs;

    //internal
    double m_max_timepoint;
    double m_min_amplitude,m_max_amplitude;
    bool m_update_scheduled;
    QImage m_image;

    //settings
    double m_hmargin,m_vmargin;

    QPointF coord2imagepix(const QPointF &p,int W,int H);
    void schedule_update();
};

MVFiringRateView::MVFiringRateView()
{
    d=new MVFiringRateViewPrivate;
    d->q=this;
    d->m_sampling_freq=30000; //Hz
    d->m_max_timepoint=0;
    d->m_min_amplitude=0;
    d->m_max_amplitude=0;
    d->m_update_scheduled=false;

    d->m_hmargin=25;
    d->m_vmargin=25;
}

MVFiringRateView::~MVFiringRateView()
{
    delete d;
}

void MVFiringRateView::setFirings(const Mda &firings)
{
    d->m_firings=firings;
}

void MVFiringRateView::setFirings(const Mda &firings)
{
    d->m_firings=firings;

    d->m_max_timepoint=compute_max(d->m_times);
    d->m_min_amplitude=compute_min(d->m_amplitudes);
    d->m_max_amplitude=compute_max(d->m_amplitudes);

    d->schedule_update();
}

void MVFiringRateView::setSamplingFreq(double ff)
{
    d->m_sampling_freq=ff;
    d->schedule_update();
}

void MVFiringRateView::setCurrentEvent(MVEvent evt)
{
    Q_UNUSED(evt)
    //finish this!
}

void MVFiringRateView::setEpochs(const QList<Epoch> &epochs)
{
    d->m_epochs=epochs;
    update();
}

void apply_kernel(int N,double *X,int kernel_rad,double *kernel) {
    double Y[N];
    for (int i=0; i<N; i++) {
        double val=0;
        for (int dd=-kernel_rad; dd<=kernel_rad; dd++) {
            int i2=i+dd;
            if ((0<=i2)&&(i2<N)) {
                val+=X[i2]*kernel[kernel_rad+dd];
            }
        }
        Y[i]=val;
    }
    for (int i=0; i<N; i++) X[i]=Y[i];
}

void smooth_grid(Mda &X,double kernel_tau) {
    int kernel_rad=(int)kernel_tau*3;
    double kernel[(kernel_rad*2+1)];
    for (int dd=-kernel_rad; dd<=kernel_rad; dd++) {
        kernel[dd+kernel_rad]=exp(-0.5*(dd*dd)/(kernel_tau*kernel_tau));
    }
    double *ptr=X.dataPtr();
    double bufx[X.N1()];
    int aa=0;
    for (int y=0; y<X.N2(); y++) {
        for (int x=0; x<X.N1(); x++) {
            bufx[x]=ptr[aa];
            aa++;
        }
        aa-=X.N1();
        apply_kernel(X.N1(),bufx,kernel_rad,kernel);
        for (int x=0; x<X.N1(); x++) {
            ptr[aa]=bufx[x];
            aa++;
        }
    }
    double bufy[X.N2()];
    for (int x=0; x<X.N1(); x++) {
        int aa=x;
        for (int y=0; y<X.N2(); y++) {
            bufy[y]=ptr[aa];
            aa+=X.N1();
        }
        apply_kernel(X.N2(),bufy,kernel_rad,kernel);
        aa=x;
        for (int y=0; y<X.N2(); y++) {
            ptr[aa]=bufy[y];
            aa+=X.N1();
        }
    }
}

void MVFiringRateView::slot_update()
{
    d->m_update_scheduled=false;

    int W=500;
    int H=(int)(this->height()*1.0/this->width()*W);

    Mda grid; grid.allocate(W,H);
    double *ptr=grid.dataPtr();
    for (int i=0; i<d->m_times.count(); i++) {
        double t0=d->m_times[i];
        double a0=d->m_amplitudes.value(i);
        QPointF pt=d->coord2imagepix(QPointF(t0,a0),W,H);
        int x=(int)pt.x();
        int y=(int)pt.y();
        if ((x>=0)&&(x<W)&&(y>=0)&&(y<H)) {
            ptr[x+W*y]++;
        }
    }

    Mda smoothed_grid=grid;
    smooth_grid(smoothed_grid,4);
    double max_density=compute_max(smoothed_grid.totalSize(),smoothed_grid.dataPtr());
    double *smoothed_ptr=smoothed_grid.dataPtr();

    d->m_image=QImage(W,H,QImage::Format_ARGB32);
    d->m_image.fill(QColor(0,0,0,0).toRgb());

    QColor white(255,255,255);
    QColor red(255,0,0);
    for (int y=0; y<H; y++) {
        for (int x=0; x<W; x++) {
            double val=ptr[x+W*y];
            if (val) {
                double pct=smoothed_ptr[x+W*y]/max_density;
                pct=sqrt(pct);
                QColor CC=get_heat_map_color(pct);
                d->m_image.setPixel(x,y,CC.rgb());
            }
        }
    }

    this->update();
    /*d->m_update_scheduled=false;
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
    d->m_view->setData(DAM,this);*/
}




void MVFiringRateView::paintEvent(QPaintEvent *evt)
{
    Q_UNUSED(evt);

    qDebug() << "paintEvent" << d->m_times.count() << d->m_amplitudes.count();

    QPainter painter(this);
    painter.fillRect(0,0,width(),height(),QColor(160,160,160));

    QRectF target=QRectF(d->m_hmargin,d->m_vmargin,width()-2*d->m_hmargin,height()-2*d->m_vmargin);

    QFont font=painter.font(); font.setPixelSize(24);
    painter.setFont(font);
    QColor epoch_color=QColor(150,150,150);
    //double margin=(d->m_max_amplitude-d->m_min_amplitude)*0.3;
    double margin=0;
    for (int i=0; i<d->m_epochs.count(); i++) {
        Epoch epoch=d->m_epochs[i];
        QPointF pt1=d->coord2imagepix(QPointF(epoch.t_begin,d->m_max_amplitude+margin),target.width(),target.height());
        QPointF pt2=d->coord2imagepix(QPointF(epoch.t_end,d->m_min_amplitude-margin),target.width(),target.height());
        QRect RR(target.x()+pt1.x(),target.y()+pt1.y(),pt2.x()-pt1.x(),pt2.y()-pt1.y());
        painter.fillRect(RR,epoch_color);
        RR.adjust(0,d->m_vmargin,0,0);
        if (fabs(d->m_max_amplitude)>fabs(d->m_min_amplitude)) {
            painter.drawText(RR,Qt::AlignTop|Qt::AlignHCenter,epoch.name);
        }
        else {
            painter.drawText(RR,Qt::AlignBottom|Qt::AlignHCenter,epoch.name);
        }
    }

    painter.drawImage(target,d->m_image.scaled(target.width(),target.height(),Qt::IgnoreAspectRatio,Qt::SmoothTransformation));
}


QPointF MVFiringRateViewPrivate::coord2imagepix(const QPointF &p,int W,int H)
{
    if ((m_max_timepoint==0)||(m_max_amplitude==m_min_amplitude)) return QPointF(0,0);
    double pctx=p.x()/m_max_timepoint;
    double pcty=(p.y()-m_min_amplitude)/(m_max_amplitude-m_min_amplitude);

    double x0=pctx*W;
    double y0=(1-pcty)*H;

    return QPointF(x0,y0);
}

void MVFiringRateViewPrivate::schedule_update()
{
    if (m_update_scheduled) return;
    m_update_scheduled=true;
    QTimer::singleShot(500,q,SLOT(slot_update()));
}
