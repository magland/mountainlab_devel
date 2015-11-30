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
#include <QList>
#include <QMessageBox>
#include <QTextBrowser>
#include "firetrackwidget.h"
#include "histogramview.h"
#include <math.h>
#include "mvcrosscorrelogramswidget.h"
#include <QSplitter>

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

	DiskArrayModel *m_clips;
	bool m_own_clips;
	bool m_unit_number;

	MVCrossCorrelogramsWidget *m_cross_correlograms_widget;
	SSTimeSeriesView *m_clips_view;

	QSplitter *m_hsplitter;
	QSplitter *m_vsplitter;

	void update_sizes();
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

	d->m_clips_view=new SSTimeSeriesView;
	d->m_clips_view->initialize();

	d->m_cross_correlograms_widget=new MVCrossCorrelogramsWidget;

	QSplitter *vsplitter=new QSplitter(Qt::Vertical);
	vsplitter->setHandleWidth(15);
	vsplitter->addWidget(d->m_clips_view);
	vsplitter->addWidget(d->m_cross_correlograms_widget);
	d->m_vsplitter=vsplitter;

	QSplitter *hsplitter=new QSplitter;
	hsplitter->setHandleWidth(15);
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

void MVUnitWidget::setClips(DiskArrayModel *C, bool own_it)
{
	if ((d->m_clips)&&(d->m_own_clips)) delete d->m_clips;
	d->m_clips=C;
	d->m_own_clips=own_it;
	d->m_clips_view->setData(C,false);
}

void MVUnitWidget::setUnitNumber(int num)
{
	d->m_unit_number=num;
}

void MVUnitWidget::updateWidgets()
{
	d->m_cross_correlograms_widget->updateWidget();
}

void MVUnitWidget::resizeEvent(QResizeEvent *evt)
{
	Q_UNUSED(evt);

	d->update_sizes();
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
}

MVUnitWidget::~MVUnitWidget()
{
	if ((d->m_raw)&&(d->m_own_raw)) delete d->m_raw;
	if ((d->m_raw_whitened)&&(d->m_own_raw_whitened)) delete d->m_raw_whitened;
	if ((d->m_clips)&&(d->m_own_clips)) delete d->m_clips;
	delete d;
}

