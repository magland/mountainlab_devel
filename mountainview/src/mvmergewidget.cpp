#include "mvmergewidget.h"
#include "diskreadmda.h"
#include "sstimeseriesview.h"
#include "correlationmatrixview.h"
#include "mvcrosscorrelogramswidget.h"

#include <QHBoxLayout>
#include <QSplitter>

class MVMergeWidgetPrivate {
public:
	MVMergeWidget *q;
	DiskReadMda m_raw;
	DiskReadMda m_clusters;
	Mda m_correlation_matrix;

	SSTimeSeriesView *m_raw_view;
	CorrelationMatrixView *m_correlation_matrix_view;
	MVCrossCorrelogramsWidget *m_cross_correlograms_widget;

	QSplitter *m_splitter1,*m_splitter2;

	Mda create_cross_correlograms_data();
	void update_sizes();
};

MVMergeWidget::MVMergeWidget(QWidget *parent) : QWidget (parent)
{
	d=new MVMergeWidgetPrivate;
	d->q=this;

	d->m_raw_view=new SSTimeSeriesView;
	d->m_raw_view->initialize();

	d->m_correlation_matrix_view=new CorrelationMatrixView;

	d->m_cross_correlograms_widget=new MVCrossCorrelogramsWidget;

	QSplitter *splitter1=new QSplitter;
	splitter1->setOrientation(Qt::Horizontal);
	d->m_splitter1=splitter1;

	QSplitter *splitter2=new QSplitter;
	splitter2->setOrientation(Qt::Vertical);
	d->m_splitter2=splitter2;

	splitter1->addWidget(d->m_correlation_matrix_view);
	splitter1->addWidget(splitter2);

	splitter2->addWidget(d->m_cross_correlograms_widget);
	splitter2->addWidget(d->m_raw_view);

	QHBoxLayout *hlayout=new QHBoxLayout;
	hlayout->addWidget(splitter1);
	this->setLayout(hlayout);
}

MVMergeWidget::~MVMergeWidget()
{
	delete d;
}

void MVMergeWidget::setRawPath(const QString &path)
{
	d->m_raw.setPath(path);
	DiskArrayModel *X=new DiskArrayModel;
	X->setPath(path);
	d->m_raw_view->setData(X,true);
}

void MVMergeWidget::setClustersPath(const QString &clusters)
{
	d->m_clusters.setPath(clusters);
	QList<double> times;
	QList<double> labels;
	for (int n=0; n<d->m_clusters.N2(); n++) {
		times << d->m_clusters.value(1,n);
		labels << d->m_clusters.value(2,n);
	}
	d->m_raw_view->setTimesLabels(times,labels);

	Mda CCD=d->create_cross_correlograms_data();
	d->m_cross_correlograms_widget->setCrossCorrelogramsData(DiskReadMda(CCD));
	d->m_cross_correlograms_widget->updateWidget();
}

void MVMergeWidget::setCorrelationMatrix(const Mda &CM)
{
	d->m_correlation_matrix=CM;
	d->m_correlation_matrix_view->setMatrix(CM);
}

void MVMergeWidget::resizeEvent(QResizeEvent *evt)
{
	Q_UNUSED(evt)
	d->update_sizes();
}

typedef QList<long> IntList;

Mda MVMergeWidgetPrivate::create_cross_correlograms_data()
{
	QList<long> times,labels;
	long L=m_clusters.N2();
	int max_dt=1500;

	printf("Setting up times and labels...\n");
	for (int n=0; n<L; n++) {
		times << (long)m_clusters.value(1,n);
		labels << (long)m_clusters.value(2,n);
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

	return ret;
}

void MVMergeWidgetPrivate::update_sizes()
{
	float W0=q->width();
	float H0=q->height();

	int W1=W0/3; if (W1<300) W1=300; if (W1>500) W1=500;
	int W2=W0-W1;
	if (W2<500) {W2=500; W1=W0-W2;}

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
