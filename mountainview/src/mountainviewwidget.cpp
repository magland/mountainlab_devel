#include "mountainviewwidget.h"
#include <QDebug>
#include <QDir>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QProgressDialog>
#include <QPushButton>
#include "sstimeserieswidget.h"
#include "sstimeseriesview.h"
#include "diskreadmda.h"
#include <QList>
#include <QTextBrowser>
#include "firetrackwidget.h"

class MountainViewWidgetPrivate {
public:
	MountainViewWidget *q;

	Mda m_primary_channels;
    Mda m_templates;
	Mda m_locations;
    DiskArrayModel *m_raw;
    Mda m_times;
    Mda m_labels;
    DiskReadMda *m_times_labels;
    int m_template_view_padding;
    QList<SSTimeSeriesView *> m_spike_template_views;
    QList<SSTimeSeriesView *> m_labeled_raw_data_views;
    QList<SSTimeSeriesView *> m_clips_views;
    int m_current_template_index;

    void connect_spike_templates_view(SSTimeSeriesView *V);
    void connect_labeled_raw_data_view(SSTimeSeriesView *V);
    void connect_clips_view(SSTimeSeriesView *V);
    void update_clips_view(SSTimeSeriesWidget *W,SSTimeSeriesView *V,int label);
};


MountainViewWidget::MountainViewWidget(QWidget *parent) : QMainWindow(parent)
{
	d=new MountainViewWidgetPrivate;
	d->q=this;

    d->m_raw=0;
    d->m_times_labels=0;
    d->m_template_view_padding=30;

    d->m_current_template_index=0;

    QVBoxLayout *VL=new QVBoxLayout;
    {
        QPushButton *B=new QPushButton("Labeled Raw Data (SpikeSpy)");
        VL->addWidget(B);
        connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_view_labeled_raw_data()));
    }
    {
        QPushButton *B=new QPushButton("Spike Templates");
        VL->addWidget(B);
        connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_view_spike_templates()));
    }
    {
        QPushButton *B=new QPushButton("Spike Clips");
        VL->addWidget(B);
        connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_view_spike_clips()));
    }
	{
		QPushButton *B=new QPushButton("Statistics");
		VL->addWidget(B);
		connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_statistics()));
	}
	{
		QPushButton *B=new QPushButton("FireTrack");
		VL->addWidget(B);
		connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_firetrack()));
	}
	{
		QPushButton *B=new QPushButton("Quit");
		VL->addWidget(B);
		connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_quit()));
	}
//    {
//        QPushButton *B=new QPushButton("Cluster View");
//        VL->addWidget(B);
//        connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_cluster_view()));
//    }
//    {
//        QPushButton *B=new QPushButton("FireTrack");
//        VL->addWidget(B);
//        connect(B,SIGNAL(clicked(bool)),this,SLOT(slot_firetrack()));
//    }


	QWidget *CW=new QWidget;
    CW->setLayout(VL);
	this->setCentralWidget(CW);

	this->setWindowFlags(Qt::WindowStaysOnTopHint);
    this->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowTitleHint);
    //this->setWindowFlags(Qt::X11BypassWindowManagerHint);
    //this->setWindowFlags(Qt::FramelessWindowHint);
}

void MountainViewWidget::setElectrodeLocations(const Mda &L)
{
	d->m_locations=L;
}

void MountainViewWidget::setTemplates(const Mda &X)
{
	qDebug() << "setTemplates" << X.N1() << X.N2() << X.N3();
	d->m_templates=X;
}

void MountainViewWidget::setPrimaryChannels(const Mda &X)
{
	d->m_primary_channels=X;
}

void MountainViewWidget::setRaw(DiskArrayModel *X)
{
    if (d->m_raw) delete d->m_raw;
    d->m_raw=X;
}

void MountainViewWidget::setTimesLabels(const Mda &times, const Mda &labels)
{
    d->m_times=times;
    d->m_labels=labels;
    int NN=d->m_times.totalSize();
    Mda TL; TL.allocate(2,NN);
    for (int ii=0; ii<NN; ii++) {
        TL.setValue(d->m_times.value1(ii),0,ii);
        TL.setValue(d->m_labels.value1(ii),1,ii);
    }

    d->m_times_labels=new DiskReadMda(TL);
}

void MountainViewWidget::resizeEvent(QResizeEvent *evt)
{
	Q_UNUSED(evt);
	{
		float W0=this->width()*0.15;
		W0=qMax(W0,80.0F);
		W0=qMin(W0,300.0F);
        //d->m_waveform_list->setFixedWidth(W0);
	}
	{
		float W0=this->width()*0.2;
		W0=qMax(W0,150.0F);
		W0=qMin(W0,350.0F);
        //d->m_plot->setFixedWidth(W0);
    }
}

void MountainViewWidget::slot_view_labeled_raw_data()
{
    if (!d->m_raw) return;
    SSTimeSeriesWidget *W=new SSTimeSeriesWidget;
    SSTimeSeriesView *V=new SSTimeSeriesView;
    W->addView(V);
    V->setData(d->m_raw,false);
    if (d->m_times_labels) {
        V->setLabels(d->m_times_labels,false);
    }
    V->initialize();

    W->show();
    W->setAttribute(Qt::WA_DeleteOnClose);
    W->move(this->topLevelWidget()->geometry().bottomLeft()+QPoint(0,50));
    W->resize(1600,700);

    d->connect_labeled_raw_data_view(V);
}

void MountainViewWidget::slot_view_spike_templates()
{
    SSTimeSeriesWidget *W=new SSTimeSeriesWidget;
    SSTimeSeriesView *V=new SSTimeSeriesView;
    W->addView(V);
    DiskArrayModel *data=new DiskArrayModel;
    int M=d->m_templates.N1();
    int T=d->m_templates.N2();
    int K=d->m_templates.N3();
    int padding=d->m_template_view_padding;
    Mda templates_formatted; templates_formatted.allocate(M,(T+padding)*K);
	qDebug() << "templates_formatted:" << templates_formatted.N1() << templates_formatted.N2() << templates_formatted.N3() << M << T << K;
    Mda TL; TL.allocate(2,K);
    for (int k=0; k<K; k++) {
        TL.setValue((T+padding)*k+T/2,0,k);
        TL.setValue(k+1,1,k);
        for (int t=0; t<T; t++) {
            for (int m=0; m<M; m++) {
                templates_formatted.setValue(d->m_templates.value(m,t,k),m,t+(T+padding)*k);
            }
        }
    }
    data->setFromMda(templates_formatted);
    V->plot()->setShowMarkerLines(false);
    V->setData(data,true);
    V->setLabels(new DiskReadMda(TL),true);
    V->initialize();

    W->show();
	W->setAttribute(Qt::WA_DeleteOnClose);
	W->move(this->topLevelWidget()->geometry().topRight()+QPoint(300,-100));
    W->resize(800,400);

    d->connect_spike_templates_view(V);
}

void MountainViewWidget::slot_view_spike_clips()
{
    if (!d->m_raw) return;

    QString txt=QInputDialog::getText(this,"View Spike Clips","Label number:",QLineEdit::Normal,"1");
    int label=txt.toInt();
    if (label<=0) return;

    SSTimeSeriesWidget *W=new SSTimeSeriesWidget;
    SSTimeSeriesView *V=new SSTimeSeriesView;
    V->setClipMode(true);
    V->setProperty("fixed_clipsize",true);
	V->setVerticalZoomFactor(0.5);
    W->addView(V);
	d->update_clips_view(W,V,label);

    W->show();
    W->setAttribute(Qt::WA_DeleteOnClose);
    W->move(this->topLevelWidget()->geometry().topRight()+QPoint(500,400));
    W->resize(800,400);

    d->connect_clips_view(V);
}

void MountainViewWidget::slot_firetrack() {
	FireTrackWidget *W=new FireTrackWidget;
	W->setElectrodeLocations(d->m_locations);
	W->setWaveforms(d->m_templates);
	W->electrodeArrayWidget()->setShowChannelNumbers(true);

	W->show();
	W->setAttribute(Qt::WA_DeleteOnClose);
	W->move(this->topLevelWidget()->geometry().topRight()+QPoint(500,400));
	W->resize(800,400);
}

void MountainViewWidget::slot_cluster_view()
{

}

struct SpikeStats {
	int count;
};

QString read_text_file(QString path) {
	QFile FF(path);
	if (!FF.open(QFile::Text|QFile::ReadOnly)) {
		return "";
	}
	QString ret=QString(FF.readAll());
	FF.close();
	return ret;
}

void MountainViewWidget::slot_statistics()
{
	QList<SpikeStats> spike_stats;
	int num=d->m_times.totalSize();
	for (int ii=0; ii<num; ii++) {
		//int t0=d->m_times.value1(ii);
		int label0=d->m_labels.value1(ii);
		while (label0>=spike_stats.count()) {
			SpikeStats X;
			X.count=0;
			spike_stats << X;
		}
		spike_stats[label0].count++;
	}

	QString statistics_css=read_text_file(":/statistics.css");

	QString html;
	html+="<html>\n";
	html+=QString("<head>\n");
	html+=QString("<style>\n");
	html+=QString("%1\n\n").arg(statistics_css);
	html+=QString("</style>\n");
	html+=QString("</head>\n");
	html+=QString("<table>\n");
	html+=QString("<tr><th>Template</th><th>Primary Channel</th><th># Spikes</th><th>Spikes per minute</th></tr>");
	for (int k=1; k<spike_stats.count(); k++) {
		SpikeStats X=spike_stats[k];
		int primary_channel=(int)d->m_primary_channels.value1(k-1);
		if (X.count!=0) {
			float frequency=X.count*1.0/(d->m_raw->size(1)*1.0/30000/60);
			html+=QString("<tr><td>%1</td><td>%2</td><td>%3</td><td>%4</td></tr>\n")
					.arg(k)
					.arg(primary_channel)
					.arg(X.count)
					.arg(frequency,0,'f',1,' ');
		}
	}
	html+=QString("</table>\n");
	html+="</html>\n";

	{
		QTextBrowser *W=new QTextBrowser;
		W->setAttribute(Qt::WA_DeleteOnClose);
		W->setHtml(html);
		W->show();

		W->resize(600,800);
		W->move(this->geometry().topRight()+QPoint(100,50));
	}
}

void MountainViewWidget::slot_quit()
{
	qApp->quit();
}

void MountainViewWidget::slot_spike_templates_x_changed()
{
    SSTimeSeriesView *V=(SSTimeSeriesView *)sender();
    int x=V->currentX();
    int T=d->m_templates.N2();
    int k=x/(T+d->m_template_view_padding) + 1;
    d->m_current_template_index=k;

    //if (d->m_clips_widget->isVisible()) {
    //    d->update_clips_view();
    //}
}

void MountainViewWidget::slot_object_destroyed(QObject *obj)
{
    d->m_spike_template_views.removeOne((SSTimeSeriesView *)obj);
    d->m_labeled_raw_data_views.removeOne((SSTimeSeriesView *)obj);
    d->m_clips_views.removeOne((SSTimeSeriesView *)obj);
}

MountainViewWidget::~MountainViewWidget()
{
    if (d->m_raw) delete d->m_raw;
    if (d->m_times_labels) delete d->m_times_labels;
	delete d;
}

void MountainViewWidgetPrivate::connect_spike_templates_view(SSTimeSeriesView *V)
{
    QObject::connect(V,SIGNAL(currentXChanged()),q,SLOT(slot_spike_templates_x_changed()));
    QObject::connect(V,SIGNAL(destroyed(QObject*)),q,SLOT(slot_object_destroyed(QObject*)));
    m_spike_template_views << V;
}

void MountainViewWidgetPrivate::connect_labeled_raw_data_view(SSTimeSeriesView *V)
{
    QObject::connect(V,SIGNAL(destroyed(QObject*)),q,SLOT(slot_object_destroyed(QObject*)));
    m_labeled_raw_data_views << V;
}

void MountainViewWidgetPrivate::connect_clips_view(SSTimeSeriesView *V)
{
    QObject::connect(V,SIGNAL(destroyed(QObject*)),q,SLOT(slot_object_destroyed(QObject*)));
    m_clips_views << V;
}

Mda extract_clips(DiskArrayModel *X,const Mda &times,const Mda &labels,int label) {
    Mda clips;

    QList<int> times0;
    for (int i=0; i<labels.totalSize(); i++) {
        if (labels.value1(i)==label) {
            times0 << (int)times.value1(i);
        }
    }

    int M=X->size(0);
    int N=X->size(1);
    int T=100;
    int NC=times0.count();

    int dtmin=-T/2;
    int dtmax=dtmin+T-1;

    clips.allocate(M,T,NC);
    for (int ii=0; ii<NC; ii++) {
        int t0=times0[ii];
        if ((t0+dtmin>=0)&&(t0+dtmax<N)) {
            for (int t=0; t<T; t++) {
                int t1=t0+dtmin+t;
                for (int m=0; m<M; m++) {
                    clips.setValue(X->value(m,t1),m,t,ii);
                }
            }
        }
    }
    return clips;
}

Mda format_clips(const Mda &clips,int padding) {
    int M=clips.N1();
    int T=clips.N2();
    int NC=clips.N3();
    Mda clips2;
    clips2.allocate(M,(T+padding)*NC);
    for (int i=0; i<NC; i++) {
        for (int t=0; t<T; t++) {
            for (int m=0; m<M; m++) {
                clips2.setValue(clips.value(m,t,i),m,t+(T+padding)*i);
            }
        }
    }
    return clips2;
}

void MountainViewWidgetPrivate::update_clips_view(SSTimeSeriesWidget *W,SSTimeSeriesView *V,int label)
{
    QProgressDialog dlg;
	dlg.setWindowTitle(QString("Extracting clips for template %1").arg(label));
    dlg.setRange(0,100);
    dlg.show();
	dlg.setLabelText(QString("Extracting clips for template %1...").arg(label));
    dlg.setValue(0); dlg.repaint(); qApp->processEvents();
    Mda clips=extract_clips(m_raw,m_times,m_labels,label);
    dlg.setLabelText("Formatting clips...");
    dlg.setValue(50); dlg.repaint(); qApp->processEvents();
    Mda clips2=format_clips(clips,m_template_view_padding);
    DiskArrayModel *MM=new DiskArrayModel;
    MM->setFromMda(clips2);
    dlg.setLabelText("Initializing...");
    dlg.setValue(75); dlg.repaint(); qApp->processEvents();
    V->setData(MM,true);
    V->initialize();
    W->setClipData(clips);
    W->setWindowTitle(QString("Spike Clips -- template %1 -- %2 spikes").arg(label).arg(clips.N3()));
}
