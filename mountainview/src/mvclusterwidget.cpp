#include "mvclusterwidget.h"
#include "mvclusterview.h"
#include <QHBoxLayout>
#include <QLabel>
#include <QList>
#include "mvclipsview.h"
#include "msutils.h"

class MVClusterWidgetPrivate {
public:
	MVClusterWidget *q;
	QList<MVClusterView *> m_views;
	MVClipsView *m_clips_view;
	QLabel *m_info_bar;
	DiskReadMda m_raw;
	int m_clip_size;
	QList<double> m_outlier_scores;

	void connect_view(MVClusterView *V);
	void update_clips_view();
	int current_event_index();
};

MVClusterWidget::MVClusterWidget()
{
	d=new MVClusterWidgetPrivate;
	d->q=this;

	d->m_clip_size=200;

	d->m_clips_view=new MVClipsView;

	{
		MVClusterView *X=new MVClusterView;
		X->setMode(MVCV_MODE_HEAT_DENSITY);
		d->m_views << X;
	}
	{
		MVClusterView *X=new MVClusterView;
		X->setMode(MVCV_MODE_LABEL_COLORS);
		d->m_views << X;
	}

	QHBoxLayout *hlayout=new QHBoxLayout;
	this->setLayout(hlayout);

	QVBoxLayout *vlayout=new QVBoxLayout;
	vlayout->addWidget(d->m_clips_view);
	d->m_info_bar=new QLabel; d->m_info_bar->setFixedHeight(20);
	vlayout->addWidget(d->m_info_bar);
	hlayout->addLayout(vlayout);
	d->m_clips_view->setFixedWidth(250);

	foreach (MVClusterView *V,d->m_views) {
		hlayout->addWidget(V);
	}

	foreach (MVClusterView *V,d->m_views) {
		d->connect_view(V);
	}
}

MVClusterWidget::~MVClusterWidget()
{
	delete d;
}

void MVClusterWidget::setData(const Mda &X)
{
	foreach (MVClusterView *V,d->m_views) {
		V->setData(X);
	}
}

void MVClusterWidget::setTimes(const QList<double> &times)
{
	foreach (MVClusterView *V,d->m_views) {
		V->setTimes(times);
	}
}

void MVClusterWidget::setLabels(const QList<int> &labels)
{
	foreach (MVClusterView *V,d->m_views) {
		V->setLabels(labels);
	}
}

void MVClusterWidget::setOutlierScores(const QList<double> &outlier_scores)
{
	d->m_outlier_scores=outlier_scores;
}

void MVClusterWidget::setCurrentEvent(const MVEvent &evt)
{
	foreach (MVClusterView *V,d->m_views) {
		V->setCurrentEvent(evt);
	}
	d->update_clips_view();
}

void MVClusterWidget::setClipSize(int clip_size)
{
	d->m_clip_size=clip_size;
}

void MVClusterWidget::setRaw(const DiskReadMda &X)
{
	d->m_raw=X;
	d->update_clips_view();
}

MVEvent MVClusterWidget::currentEvent()
{
	return d->m_views[0]->currentEvent();
}

void MVClusterWidget::slot_view_current_event_changed()
{
	MVClusterView *V0=(MVClusterView *)sender();
	this->setCurrentEvent(V0->currentEvent());
	emit currentEventChanged();
}

void MVClusterWidget::slot_view_transformation_changed()
{
	MVClusterView *V0=(MVClusterView *)sender();
	AffineTransformation T=V0->transformation();
	foreach (MVClusterView *V,d->m_views) {
		V->setTransformation(T);
	}
}


void MVClusterWidgetPrivate::connect_view(MVClusterView *V)
{
	QObject::connect(V,SIGNAL(currentEventChanged()),q,SLOT(slot_view_current_event_changed()));
	QObject::connect(V,SIGNAL(transformationChanged()),q,SLOT(slot_view_transformation_changed()));
}

void MVClusterWidgetPrivate::update_clips_view()
{
	MVEvent evt=q->currentEvent();
	QString info_txt;
	if (evt.time>=0) {
		QList<long> times; times << (long)evt.time;
		Mda clip0=extract_clips(m_raw,times,m_clip_size);
		double ppp=m_outlier_scores.value(current_event_index());
		if (ppp) {
			info_txt=QString("Outlier score: %1").arg(ppp);
		}
		m_clips_view->setClips(clip0);
	}
	else {
		m_clips_view->setClips(Mda());
	}
	m_info_bar->setText(info_txt);
}

int MVClusterWidgetPrivate::current_event_index()
{
	return m_views[0]->currentEventIndex();
}
