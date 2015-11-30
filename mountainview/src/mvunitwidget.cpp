#include "mvunitwidget.h"
#include <QDebug>
#include <QDir>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QProgressDialog>
#include <QPushButton>
#include "sstimeserieswidget.h"
#include "sstimeseriesview.h"
#include "diskreadmda.h"
#include <QLabel>
#include <QList>
#include <QMessageBox>
#include <QTextBrowser>
#include <QTime>
#include <QTimer>
#include "firetrackwidget.h"
#include "histogramview.h"
#include <math.h>
#include "mvcrosscorrelogramswidget.h"
#include <QSplitter>
#include "cvwidget.h"
#include "get_principal_components.h"

class MVUnitWidgetPrivate {
public:
	MVUnitWidget *q;

	Mda m_primary_channels;
	Mda m_templates,m_templates_whitened;
	Mda m_locations;
	DiskArrayModel *m_raw,*m_raw_whitened;
	bool m_own_raw,m_own_raw_whitened;
	Mda m_times;
	Mda m_labels;
	QString m_cross_correlograms_path;
	int m_current_clip_number;

	DiskArrayModel *m_clips;
	bool m_own_clips;
	bool m_unit_number;

	MVCrossCorrelogramsWidget *m_cross_correlograms_widget;
	SSTimeSeriesView *m_clips_view;
	SSTimeSeriesView *m_template_view;
	CVWidget *m_cluster_widget;

	QSplitter *m_hsplitter;
	QSplitter *m_vsplitter;
	QSplitter *m_left_vsplitter;
	QLabel *m_status_label;

	void update_sizes();
	void set_status_text(const QString &txt);
	void set_current_clip_number(int num);
};


MVUnitWidget::MVUnitWidget(QWidget *parent) : QMainWindow(parent)
{
	d=new MVUnitWidgetPrivate;
	d->q=this;

	d->m_raw=0;
	d->m_raw_whitened=0;
	d->m_own_raw=d->m_own_raw_whitened=false;

	d->m_clips=0;
	d->m_own_clips=false;
	d->m_unit_number=0;
	d->m_current_clip_number=-1;

	d->m_clips_view=new SSTimeSeriesView;
	d->m_clips_view->initialize();
	connect(d->m_clips_view,SIGNAL(currentXChanged()),this,SLOT(slot_clips_view_current_x_changed()));

	d->m_template_view=new SSTimeSeriesView;
	d->m_template_view->initialize();

	d->m_cluster_widget=new CVWidget;
	d->m_cluster_widget->setNumDataPointsToSelect(1);
	connect(d->m_cluster_widget,SIGNAL(selectedDataPointsChanged()),this,SLOT(slot_selected_data_points_changed()));

	d->m_cross_correlograms_widget=new MVCrossCorrelogramsWidget;

	d->m_status_label=new QLabel;

	{
		QSplitter *vsplitter=new QSplitter(Qt::Vertical);
		vsplitter->setHandleWidth(15);
		vsplitter->addWidget(d->m_clips_view);
		vsplitter->addWidget(d->m_cross_correlograms_widget);
		d->m_vsplitter=vsplitter;
	}

	{
		QSplitter *vsplitter=new QSplitter(Qt::Vertical);
		vsplitter->setHandleWidth(15);
		vsplitter->addWidget(d->m_template_view);
		vsplitter->addWidget(d->m_cluster_widget);
		d->m_left_vsplitter=vsplitter;
	}

	QWidget *left_panel=new QWidget;
	QVBoxLayout *left_panel_layout=new QVBoxLayout;
	left_panel->setLayout(left_panel_layout);
	left_panel_layout->addWidget(d->m_left_vsplitter);
	left_panel_layout->addWidget(d->m_status_label);

	QSplitter *hsplitter=new QSplitter;
	hsplitter->setHandleWidth(15);
	hsplitter->addWidget(left_panel);
	hsplitter->addWidget(d->m_vsplitter);
	d->m_hsplitter=hsplitter;

	QWidget *CW=new QWidget;
	QHBoxLayout *hlayout=new QHBoxLayout;
	hlayout->addWidget(hsplitter);
	CW->setLayout(hlayout);
	this->setCentralWidget(CW);

	//this->setWindowFlags(Qt::Window|Qt::WindowStaysOnTopHint);
	//this->setWindowFlags(Qt::Window|Qt::Tool);
	//this->setWindowFlags(this->windowFlags()|Qt::CustomizeWindowHint);
	//this->setWindowFlags(this->windowFlags() ^ Qt::WindowCloseButtonHint);
}

void MVUnitWidget::setElectrodeLocations(const Mda &L)
{
	d->m_locations=L;
}

void MVUnitWidget::setTemplates(const Mda &X)
{
	d->m_templates=X;
}
void MVUnitWidget::setTemplatesWhitened(const Mda &X)
{
	d->m_templates_whitened=X;
}

void MVUnitWidget::setPrimaryChannels(const Mda &X)
{
	d->m_primary_channels=X;
}

void MVUnitWidget::setRaw(DiskArrayModel *X,bool own_it)
{
	if ((d->m_raw)&&(d->m_own_raw)) delete d->m_raw;
	d->m_raw=X;
	d->m_own_raw=own_it;
}

void MVUnitWidget::setRawWhitened(DiskArrayModel *X,bool own_it)
{
	if ((d->m_raw_whitened)&&(d->m_own_raw_whitened)) delete d->m_raw_whitened;
	d->m_raw_whitened=X;
	d->m_own_raw_whitened=own_it;
}

void MVUnitWidget::setTimesLabels(const Mda &times, const Mda &labels)
{
	d->m_times=times;
	d->m_labels=labels;
//	int NN=d->m_times.totalSize();
//	Mda TL; TL.allocate(2,NN);
//	for (int ii=0; ii<NN; ii++) {
//		TL.setValue(d->m_times.value1(ii),0,ii);
//		TL.setValue(d->m_labels.value1(ii),1,ii);
//	}

//	d->m_clips_view->setLabels(new DiskReadMda(TL),true);
}

void MVUnitWidget::setCrossCorrelogramsPath(const QString &path)
{
	d->m_cross_correlograms_path=path;
	d->m_cross_correlograms_widget->setCrossCorrelogramsPath(path);
}

Mda compute_mean_waveform(DiskArrayModel *C) {
	Mda ret;
	if (!C->dim3()) return ret;
	int M=C->size(0);
	int T=C->size(1)/C->dim3();
	int NC=C->dim3();
	if (!NC) return ret;

	double sum[M*T];
	for (int ii=0; ii<M*T; ii++) sum[ii]=0;
	for (int c=0; c<NC; c++) {
		if ((c%100==0)||(c==NC-1)) {
			qApp->processEvents();
			//int pct=(int)(c*1.0/NC*100);
			//printf("Computing mean waveform...%d/%d (%d%%)\n",c,NC,pct);
		}
		int ii=0;
		Mda tmp=C->loadData(1,T*c,T*(c+1));
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				sum[ii]+=tmp.value(m,t);
				ii++;
			}
		}
	}
	ret.allocate(M,T);
	{
		int ii=0;
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				ret.setValue(sum[ii]/NC,m,t);
				ii++;
			}
		}
	}
	return ret;
}

Mda compute_features(DiskArrayModel *C) {
	Mda ret;
	if (!C->dim3()) return ret;
	int M=C->size(0);
	int T=C->size(1)/C->dim3();
	int NC=C->dim3();
	if (!NC) return ret;

	Mda X=C->loadData(1,0,T*NC);
	ret.allocate(3,NC);
	get_pca_features(M*T,NC,3,ret.dataPtr(),X.dataPtr());

	return ret;
}

void MVUnitWidget::setClips(DiskArrayModel *C, bool own_it)
{
	if ((d->m_clips)&&(d->m_own_clips)) delete d->m_clips;
	d->m_clips=C;
	d->m_own_clips=own_it;
	d->m_clips_view->setData(C,false);
	int max0=10000;
	if (C->size(1)>max0) {
		d->m_clips_view->setXRange(vec2(0,max0-1));
	}
}

void MVUnitWidget::setUnitNumber(int num)
{
	d->m_unit_number=num;
	d->m_cross_correlograms_widget->setBaseUnit(num);
}

void MVUnitWidget::updateWidgets()
{
	d->m_cross_correlograms_widget->updateWidget();
	QTimer::singleShot(100,this,SLOT(slot_compute_template()));
}

void MVUnitWidget::resizeEvent(QResizeEvent *evt)
{
	Q_UNUSED(evt);

	d->update_sizes();
}

void MVUnitWidget::slot_compute_template()
{
	d->set_status_text("Computing Template...");
	Mda template0=compute_mean_waveform(d->m_clips);
	d->set_status_text("Ready.");
	DiskArrayModel *template_data=new DiskArrayModel;
	template_data->setFromMda(template0);
	d->m_template_view->setData(template_data,true);

	d->set_status_text("Computing Features...");
	Mda features0=compute_features(d->m_clips);
	d->set_status_text("Ready.");
	d->m_cluster_widget->setFeatures(features0);
	d->m_cluster_widget->autoSetRange();
	d->m_cluster_widget->refresh();
}

void MVUnitWidget::slot_clips_view_current_x_changed()
{
	if (!d->m_clips->dim3()) return;
	//int M=d->m_clips->size(0);
	int T=d->m_clips->size(1)/d->m_clips->dim3();
	int NC=d->m_clips->dim3();
	if (!NC) return;

	int x=d->m_clips_view->currentX();
	int clip_number=x/T;
	d->set_current_clip_number(clip_number);
}

void MVUnitWidget::slot_selected_data_points_changed()
{
	QList<int> L=d->m_cluster_widget->selectedDataPointIndices();
	if (L.count()!=1) return;
	d->set_current_clip_number(L[0]);
}



void MVUnitWidgetPrivate::update_sizes()
{
	float W0=q->width();
	float H0=q->height();

	int W1=W0/3; if (W1<300) W1=300; if (W1>500) W1=500;
	int W2=W0-W1;
	if (W2<500) {W2=500; W1=W0-W2;}

	int H1=H0/3;
	int H2=H0/3;
	int H3=H0-H1-H2;

	{
		QList<int> sizes; sizes << W1 << W2;
		m_hsplitter->setSizes(sizes);
	}
	{
		QList<int> sizes; sizes << H1 << H2 << H3;
		m_vsplitter->setSizes(sizes);
	}

	int H1_left=H0/2;
	int H2_left=H0/2;
	{
		QList<int> sizes; sizes << H1_left << H2_left;
		m_left_vsplitter->setSizes(sizes);
	}

}

void MVUnitWidgetPrivate::set_status_text(const QString &txt)
{
	m_status_label->setText(txt);
}

void MVUnitWidgetPrivate::set_current_clip_number(int num)
{
	if (m_current_clip_number==num) return;
	m_current_clip_number=num;
	QList<int> L; L << m_current_clip_number;
	qDebug() << "setSelectedDataPointIndices" << L;
	m_cluster_widget->setSelectedDataPointIndices(L);

	{
		if (!m_clips->dim3()) return;
		//int M=d->m_clips->size(0);
		int T=m_clips->size(1)/m_clips->dim3();
		int NC=m_clips->dim3();
		if (!NC) return;

		int x=m_clips_view->currentX();
		int clip_number2=x/T;
		if (clip_number2!=num) {
			m_clips_view->setCurrentX(num*T+T/2);
		}
	}
}

MVUnitWidget::~MVUnitWidget()
{
	if ((d->m_raw)&&(d->m_own_raw)) delete d->m_raw;
	if ((d->m_raw_whitened)&&(d->m_own_raw_whitened)) delete d->m_raw_whitened;
	if ((d->m_clips)&&(d->m_own_clips)) delete d->m_clips;
	delete d;
}
