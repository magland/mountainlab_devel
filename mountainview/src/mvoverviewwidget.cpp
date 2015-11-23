#include "mvoverviewwidget.h"
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
#include <QMessageBox>
#include <QTextBrowser>
#include "firetrackwidget.h"
#include "histogramview.h"
#include <math.h>
#include "mvstatisticswidget.h"
#include "mvcrosscorrelogramswidget.h"
#include <QSplitter>

class MVOverviewWidgetPrivate {
public:
	MVOverviewWidget *q;

	Mda m_primary_channels;
	Mda m_templates,m_templates_whitened;
	Mda m_locations;
	DiskArrayModel *m_raw,*m_raw_whitened;
	Mda m_times;
	Mda m_labels;
	DiskReadMda *m_times_labels;
	int m_template_view_padding;
	QString m_cross_correlograms_path;
	int m_current_unit_number;

	SSTimeSeriesView *m_spike_templates_view;
	MVStatisticsWidget *m_statistics_widget;
	MVCrossCorrelogramsWidget *m_cross_correlograms_widget;
	SSTimeSeriesView *m_labeled_raw_view;

	QSplitter *m_hsplitter;
	QSplitter *m_vsplitter;

	void update_spike_templates();
	void update_sizes();
	void set_current_unit(int num);

};


MVOverviewWidget::MVOverviewWidget(QWidget *parent) : QMainWindow(parent)
{
	d=new MVOverviewWidgetPrivate;
	d->q=this;

	d->m_raw=0;
	d->m_raw_whitened=0;
	d->m_times_labels=0;
	d->m_template_view_padding=30;

	d->m_spike_templates_view=new SSTimeSeriesView;
	d->m_spike_templates_view->initialize();
	connect(d->m_spike_templates_view,SIGNAL(currentXChanged()),this,SLOT(slot_spike_templates_current_x_changed()));

	d->m_labeled_raw_view=new SSTimeSeriesView;
	d->m_labeled_raw_view->initialize();

	d->m_statistics_widget=new MVStatisticsWidget;
	connect(d->m_statistics_widget,SIGNAL(currentUnitChanged()),this,SLOT(slot_statistics_widget_current_unit_changed()));

	d->m_cross_correlograms_widget=new MVCrossCorrelogramsWidget;
	connect(d->m_cross_correlograms_widget,SIGNAL(currentUnitChanged()),this,SLOT(slot_cross_correlograms_current_unit_changed()));

	QSplitter *vsplitter=new QSplitter(Qt::Vertical);
	vsplitter->setHandleWidth(15);
	vsplitter->addWidget(d->m_spike_templates_view);
	vsplitter->addWidget(d->m_cross_correlograms_widget);
	vsplitter->addWidget(d->m_labeled_raw_view);
	d->m_vsplitter=vsplitter;

	QSplitter *hsplitter=new QSplitter;
	hsplitter->setHandleWidth(15);
	hsplitter->addWidget(d->m_statistics_widget);
	hsplitter->addWidget(vsplitter);
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

void MVOverviewWidget::setElectrodeLocations(const Mda &L)
{
	d->m_locations=L;
}

void MVOverviewWidget::setTemplates(const Mda &X)
{
	d->m_templates=X;
}
void MVOverviewWidget::setTemplatesWhitened(const Mda &X)
{
	d->m_templates_whitened=X;
}

void MVOverviewWidget::setPrimaryChannels(const Mda &X)
{
	d->m_primary_channels=X;
	d->m_statistics_widget->setPrimaryChannels(X);
}

void MVOverviewWidget::setRaw(DiskArrayModel *X)
{
	if (d->m_raw) delete d->m_raw;
	d->m_raw=X;
	d->m_statistics_widget->setRaw(X);
	d->m_labeled_raw_view->setData(X);
}

void MVOverviewWidget::setRawWhitened(DiskArrayModel *X)
{
	if (d->m_raw_whitened) delete d->m_raw_whitened;
	d->m_raw_whitened=X;
}

void MVOverviewWidget::setTimesLabels(const Mda &times, const Mda &labels)
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
	d->m_statistics_widget->setTimesLabels(times,labels);
	d->m_labeled_raw_view->setLabels(new DiskReadMda(TL),true);
}

void MVOverviewWidget::setCrossCorrelogramsPath(const QString &path)
{
	d->m_cross_correlograms_path=path;
	d->m_cross_correlograms_widget->setCrossCorrelogramsPath(path);
}

void MVOverviewWidget::updateWidgets()
{
	d->update_spike_templates();
	d->m_statistics_widget->updateStatistics();
	d->m_cross_correlograms_widget->updateWidget();
}

void MVOverviewWidget::resizeEvent(QResizeEvent *evt)
{
	Q_UNUSED(evt);

	d->update_sizes();
}

void MVOverviewWidget::slot_spike_templates_current_x_changed()
{
	int x=d->m_spike_templates_view->currentX();
	int unit_number=x/(d->m_templates.N2()+d->m_template_view_padding)+1;
	d->set_current_unit(unit_number);
}

void MVOverviewWidget::slot_cross_correlograms_current_unit_changed()
{
	int num=d->m_cross_correlograms_widget->currentUnit();
	d->set_current_unit(num);
}

void MVOverviewWidget::slot_statistics_widget_current_unit_changed()
{
	int num=d->m_statistics_widget->currentUnit();
	d->set_current_unit(num);
}

void MVOverviewWidgetPrivate::update_spike_templates()
{
	SSTimeSeriesView *V=m_spike_templates_view;
	DiskArrayModel *data=new DiskArrayModel;
	int M=m_templates.N1();
	int T=m_templates.N2();
	int K=m_templates.N3();
	int padding=m_template_view_padding;
	Mda templates_formatted; templates_formatted.allocate(M,(T+padding)*K);
	Mda TL; TL.allocate(2,K);
	Mda *templates=&m_templates;
	//if (sender()->property("whitened").toBool()) {
	//	templates=&m_templates_whitened;
	//}
	//QList<int> order=get_template_sort_order(*templates);
	for (int k=0; k<K; k++) {
		TL.setValue((T+padding)*k+T/2,0,k);
		//TL.setValue(order[k]+1,1,k);
		TL.setValue(k+1,1,k);
		for (int t=0; t<T; t++) {
			for (int m=0; m<M; m++) {
				//templates_formatted.setValue(templates->value(m,t,order[k]),m,t+(T+padding)*k);
				templates_formatted.setValue(templates->value(m,t,k),m,t+(T+padding)*k);
			}
		}
	}
	data->setFromMda(templates_formatted);
	V->plot()->setShowMarkerLines(false);
	V->setData(data,true);
	V->setVerticalZoomFactor(0.6);
	V->setLabels(new DiskReadMda(TL),true);
	V->initialize();
}

void MVOverviewWidgetPrivate::update_sizes()
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
}

void MVOverviewWidgetPrivate::set_current_unit(int num)
{
	m_current_unit_number=num;
	m_cross_correlograms_widget->setCurrentUnit(num);
	m_statistics_widget->setCurrentUnit(num);

	int x=m_spike_templates_view->currentX();
	int factor=(m_templates.N2()+m_template_view_padding);
	int unit_number=x/factor+1;
	if (unit_number!=num) {
		m_spike_templates_view->setCurrentX(factor*(num-1)+m_templates.N2()/2);
	}


}



MVOverviewWidget::~MVOverviewWidget()
{
	if (d->m_raw) delete d->m_raw;
	if (d->m_raw_whitened) delete d->m_raw_whitened;
	if (d->m_times_labels) delete d->m_times_labels;
	delete d;
}

