#include "mvoverview2widget.h"
#include "diskreadmda.h"
#include "sstimeseriesview.h"
#include "mvcrosscorrelogramswidget.h"
#include "mvoverview2widgetcontrolpanel.h"

#include <QHBoxLayout>
#include <QSplitter>
#include <QTabWidget>

class MVOverview2WidgetPrivate {
public:
	MVOverview2Widget *q;
	DiskReadMda m_raw;
	DiskReadMda m_firings_original;
	Mda m_firings;

	SSTimeSeriesView *m_raw_view;
	MVCrossCorrelogramsWidget *m_cross_correlograms_widget;
	SSTimeSeriesView *m_templates_widget;
	MVOverview2WidgetControlPanel *m_control_panel;

	QSplitter *m_splitter1,*m_splitter2;

	Mda create_cross_correlograms_data();
	Mda create_templates_data();
	void update_sizes();
	void update_cross_correlograms();
	void update_templates();
	void do_amplitude_split();
};

MVOverview2Widget::MVOverview2Widget(QWidget *parent) : QWidget (parent)
{
	d=new MVOverview2WidgetPrivate;
	d->q=this;

	d->m_raw_view=new SSTimeSeriesView;
	d->m_raw_view->initialize();

	d->m_cross_correlograms_widget=new MVCrossCorrelogramsWidget;

	d->m_control_panel=new MVOverview2WidgetControlPanel;
	connect(d->m_control_panel,SIGNAL(signalButtonClicked(QString)),this,SLOT(slot_control_panel_button_clicked(QString)));

	d->m_templates_widget=new SSTimeSeriesView;
	d->m_templates_widget->initialize();

	QSplitter *splitter1=new QSplitter;
	splitter1->setOrientation(Qt::Horizontal);
	d->m_splitter1=splitter1;

	QSplitter *splitter2=new QSplitter;
	splitter2->setOrientation(Qt::Vertical);
	d->m_splitter2=splitter2;

	splitter1->addWidget(d->m_control_panel);
	splitter1->addWidget(splitter2);


	QTabWidget *TW=new QTabWidget;
	TW->addTab(d->m_cross_correlograms_widget,"Cross-Correlograms");
	TW->addTab(d->m_templates_widget,"Templates");

	splitter2->addWidget(TW);
	splitter2->addWidget(d->m_raw_view);

	QHBoxLayout *hlayout=new QHBoxLayout;
	hlayout->addWidget(splitter1);
	this->setLayout(hlayout);
}

MVOverview2Widget::~MVOverview2Widget()
{
	delete d;
}

void MVOverview2Widget::setRawPath(const QString &path)
{
	d->m_raw.setPath(path);
	DiskArrayModel *X=new DiskArrayModel;
	X->setPath(path);
	d->m_raw_view->setData(X,true);
}

void MVOverview2Widget::setFiringsPath(const QString &firings)
{
	d->m_firings_original.setPath(firings);
	d->m_firings.allocate(d->m_firings_original.N1(),d->m_firings_original.N2());
	for (int i2=0; i2<d->m_firings.N2(); i2++) {
		for (int i1=0; i1<d->m_firings.N1(); i1++) {
			d->m_firings.setValue(d->m_firings_original.value(i1,i2),i1,i2);
		}
	}
	QList<double> times;
	QList<double> labels;
	for (int n=0; n<d->m_firings.N2(); n++) {
		times << d->m_firings.value(1,n);
		labels << d->m_firings.value(2,n);
	}
	d->m_raw_view->setTimesLabels(times,labels);
}

void MVOverview2Widget::updateWidgets()
{
	d->update_cross_correlograms();
	d->update_templates();
}

void MVOverview2Widget::resizeEvent(QResizeEvent *evt)
{
	Q_UNUSED(evt)
	d->update_sizes();
}

void MVOverview2Widget::slot_control_panel_button_clicked(QString str)
{
	qDebug() << "slot_control_panel_button_clicked" << str;
	if (str=="update_cross_correlograms") {
		d->update_cross_correlograms();
	}
	else if (str=="update_templates") {
		d->update_templates();
	}
	else if (str=="amplitude_split") {
		d->do_amplitude_split();
		d->update_cross_correlograms();
		d->update_templates();
	}
}

typedef QList<long> IntList;

Mda MVOverview2WidgetPrivate::create_cross_correlograms_data()
{
	QList<long> times,labels;
	long L=m_firings.N2();
	float samplefreq=30000;
	int max_dt=(int)(m_control_panel->getParameterValue("max_dt").toInt()*samplefreq/1000);
	qDebug() << "Using max_dt = " << max_dt;

	printf("Setting up times and labels...\n");
	for (int n=0; n<L; n++) {
		times << (long)m_firings.value(1,n);
		labels << (long)m_firings.value(2,n);
	}
	int K=0;
	for (int n=0; n<labels.count(); n++) {
		if (labels[n]>K) K=labels[n];
	}

	printf("Initializing output...\n");
	QList<IntList> out;
	IntList empty_list;
	for (int k1=1; k1<=K; k1++) {
		for (int k2=1; k2<=K; k2++) {
			out << empty_list;
		}
	}

	printf("Setting time differences...\n");
	int i1=0;
	for (int i2=0; i2<L; i2++) {
		while ((i1+1<L)&&(times[i1]<times[i2]-max_dt)) i1++;
		int k2=labels[i2];
		int t2=times[i2];
		if (k2>=1) {
			for (int jj=i1; jj<i2; jj++) {
				int k1=labels[jj];
				int t1=times[jj];
				if (k1>=1) {
					out[(k1-1)+K*(k2-1)] << t2-t1;
					out[(k2-1)+K*(k1-1)] << t1-t2;
				}
			}
		}
	}

	printf("Counting...\n");
	int ct=0;
	for (int k1=1; k1<=K; k1++) {
		for (int k2=1; k2<=K; k2++) {
			ct+=out[(k1-1)+K*(k2-1)].count();
		}
	}

	printf("Creating mda...\n");
	Mda ret; ret.allocate(3,ct);

	ct=0;
	for (int k1=1; k1<=K; k1++) {
		for (int k2=1; k2<=K; k2++) {
			IntList *tmp=&out[(k1-1)+K*(k2-1)];
			for (int jj=0; jj<tmp->count(); jj++) {
				ret.setValue(k1,0,ct);
				ret.setValue(k2,1,ct);
				ret.setValue((*tmp)[jj],2,ct);
				ct++;
			}
		}
	}
	printf(".\n");

	return ret;
}

Mda MVOverview2WidgetPrivate::create_templates_data()
{
	QList<long> times,labels;
	long L=m_firings.N2();
	int M=m_raw.N1();
	int T=m_control_panel->getParameterValue("clip_size").toInt();

	printf("Setting up times and labels...\n");
	for (int n=0; n<L; n++) {
		times << (long)m_firings.value(1,n);
		labels << (long)m_firings.value(2,n);
	}
	int K=0;
	for (int n=0; n<labels.count(); n++) {
		if (labels[n]>K) K=labels[n];
	}

	int dt1=-T/2;
	int dt2=dt1+T-1;

	printf("Creating mda...\n");
	Mda ret; ret.allocate(M,T,K);
	double *retptr=ret.dataPtr();
	QList<long> counts;
	for (int k=0; k<=K; k++) counts << k;
	for (int ii=0; ii<times.count(); ii++) {
		int t=times[ii];
		int k=labels[ii];
		if (k>0) {
			if ((t-1+dt1>=0)&&(t-1+dt2<m_raw.N2())) {
				counts[k]++;
				long iii=M*T*(k-1);
				long jjj=(t-1+dt1)*M;
				for (int dt=dt1; dt<=dt2; dt++) {
					for (int m=0; m<M; m++) {
						double val=m_raw.value1(jjj);
						retptr[iii]+=val;
						iii++;
						jjj++;
					}
				}
			}
		}
	}
	for (int k=1; k<=K; k++) {
		if (counts[k]>0) {
			for (int t=0; t<T; t++) {
				for (int m=0; m<M; m++) {
					double val=ret.value(m,t,k-1)/counts[k];
					ret.setValue(val,m,t,k-1);
				}
			}
		}
	}

	qDebug() << "Created templates data:" << ret.N1() << ret.N2() << ret.N3();
	return ret;
}

void MVOverview2WidgetPrivate::update_sizes()
{
	float W0=q->width();
	float H0=q->height();

	int W1=W0/3; if (W1<150) W1=150; if (W1>400) W1=400;
	int W2=W0-W1;

	int H1=H0/2;
	int H2=H0/2;
	//int H3=H0-H1-H2;

	{
		QList<int> sizes; sizes << W1 << W2;
		m_splitter1->setSizes(sizes);
	}
	{
		QList<int> sizes; sizes << H1 << H2;
		m_splitter2->setSizes(sizes);
	}

}

void MVOverview2WidgetPrivate::update_cross_correlograms()
{
	Mda CCD=create_cross_correlograms_data();
	printf("Setting cross correlograms data...\n");
	m_cross_correlograms_widget->setCrossCorrelogramsData(DiskReadMda(CCD));
	printf("Updating CC widget...\n");
	m_cross_correlograms_widget->updateWidget();
	printf(".\n");
}

void MVOverview2WidgetPrivate::update_templates()
{
	Mda TD=create_templates_data();
	printf("Setting templates data...\n");
	DiskArrayModel *MM=new DiskArrayModel;
	MM->setFromMda(TD);
	m_templates_widget->setData(MM,true);
	printf(".\n");
}

void MVOverview2WidgetPrivate::do_amplitude_split()
{
	float shell_width=m_control_panel->getParameterValue("shell_width").toFloat();
	int min_per_shell=m_control_panel->getParameterValue("min_per_shell").toInt();
	QList<long> times;
	QList<long> labels;
	QList<double> peaks;
	for (int n=0; n<m_firings.N2(); n++) {
		times << (long)m_firings_original.value(1,n);
		labels << (long)m_firings_original.value(2,n);
		peaks << m_firings_original.value(3,n);
	}
	int K=0;
	for (int n=0; n<times.count(); n++) {
		if (labels[n]>K) K=labels[n];
	}
	QList<int> nums;
	QList<float> mins;
	QList<float> maxs;
	int MAXAMP=100;
	int BINS_SIZE=(int)(2*MAXAMP/shell_width);
	for (int k=1; k<=K; k++) {
		long bins[BINS_SIZE];
		for (int ii=0; ii<BINS_SIZE; ii++) bins[ii]=0;
		long tot_count=0;
		for (int n=0; n<times.count(); n++) {
			if (labels[n]==k) {
				int ind=(int)((peaks[n]-(-MAXAMP))/shell_width);
				if ((0<=ind)&&(ind<BINS_SIZE)) {
					bins[ind]++;
					tot_count++;
				}
			}
		}
		float tmp_min=-MAXAMP;
		int running_count=0;
		int tot_used=0;
		for (int ii=0; ii<BINS_SIZE; ii++) {
			running_count+=bins[ii];
			if ((running_count>=min_per_shell)&&(tot_count-tot_used-running_count>=min_per_shell)) {
				nums << k;
				mins << tmp_min;
				maxs << -MAXAMP+(ii+1)*shell_width;
				tmp_min=-MAXAMP+(ii+1)*shell_width;
				running_count=0;
			}
		}
		if (running_count>0) {
			nums << k;
			mins << tmp_min;
			maxs << MAXAMP;
		}
	}
	int KK=nums.count();
	m_firings.allocate(m_firings_original.N1(),m_firings_original.N2());
	for (int i2=0; i2<m_firings.N2(); i2++) {
		for (int i1=0; i1<m_firings.N1(); i1++) {
			if (i1!=2) { //don't set the labels!
				m_firings.setValue(m_firings_original.value(i1,i2),i1,i2);
			}
		}
	}

	for (int kk=0; kk<KK; kk++) {
		int k=nums[kk];
		float min0=mins[kk];
		float max0=maxs[kk];
		for (int n=0; n<times.count(); n++) {
			if (labels[n]==k) {
				if ((min0<=peaks[n])&&(peaks[n]<max0)) {
					m_firings.setValue(kk,2,n);
				}
			}
		}
	}

}
