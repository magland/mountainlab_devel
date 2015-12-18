#include "mvcrosscorrelogramswidget.h"
#include <QApplication>
#include <QGridLayout>
#include <QProgressDialog>
#include "diskreadmda.h"
#include <math.h>
#include "histogramview.h"
#include <QDebug>
#include <QKeyEvent>

class MVCrossCorrelogramsWidgetPrivate {
public:
	MVCrossCorrelogramsWidget *q;
	QString m_path;
	int m_base_unit_num;
	int m_current_unit_num;
	QList<HistogramView *> m_histogram_views;
	int m_num_columns;

	void do_highlighting();
};

MVCrossCorrelogramsWidget::MVCrossCorrelogramsWidget()
{
	d=new MVCrossCorrelogramsWidgetPrivate;
	d->q=this;
	d->m_current_unit_num=0;
	d->m_base_unit_num=0;
	d->m_num_columns=-1;

	this->setFocusPolicy(Qt::StrongFocus);
}

MVCrossCorrelogramsWidget::~MVCrossCorrelogramsWidget()
{
	delete d;
}

void MVCrossCorrelogramsWidget::setCrossCorrelogramsPath(const QString &path)
{
	d->m_path=path;
}

typedef QList<float> FloatList;

QList<FloatList> get_cross_correlogram_datas_2(DiskReadMda &X,int k) {
	int K=0;
	for (int i=0; i<X.N2(); i++) {
		int k1=(int)X.value(0,i);
		int k2=(int)X.value(1,i);
		if (k1>K) K=k1;
		if (k2>K) K=k2;
	}

	QList<FloatList> ret;
	for (int i=0; i<=K; i++) {
		ret << FloatList();
	}
	for (int i=0; i<X.N2(); i++) {
		int k1=(int)X.value(0,i);
		int k2=(int)X.value(1,i);
		if ((k1==k)||((k==0)&&(k1==k2))) {
			if ((1<=k2)&&(k2<=K)) {
				ret[k2] << X.value(2,i);
			}
		}
	}
	return ret;
}

void MVCrossCorrelogramsWidget::updateWidget()
{
	if (d->m_path.isEmpty()) return;

	int k0=d->m_base_unit_num;

	DiskReadMda X;
	X.setPath(d->m_path);

	QProgressDialog dlg;
	dlg.show();
	dlg.setLabelText("Loading cross correlograms...");
	dlg.repaint(); qApp->processEvents();
	QList<FloatList> data0=get_cross_correlogram_datas_2(X,k0);

	int K=data0.count()-1;
	int num_rows=(int)sqrt(K); if (num_rows<1) num_rows=1;
	int num_cols=(K+num_rows-1)/num_rows;
	d->m_num_columns=num_cols;

	QWidget *W=this;
	W->setAttribute(Qt::WA_DeleteOnClose);
	QGridLayout *GL=new QGridLayout;
	GL->setHorizontalSpacing(20); GL->setVerticalSpacing(0);
	GL->setMargin(0);
	W->setLayout(GL);

	for (int k1=1; k1<=K; k1++) {
		HistogramView *HV=new HistogramView;
		HV->setData(data0[k1]);
		HV->autoSetBins(50);
		int k2=k1; if (k0>=1) k2=k0;
		QString title0=QString("%1/%2").arg(k1).arg(k2);
		HV->setTitle(title0);
		GL->addWidget(HV,(k1-1)/num_cols,(k1-1)%num_cols);
		HV->setProperty("unit_number",k1);
		connect(HV,SIGNAL(clicked()),this,SLOT(slot_histogram_view_clicked()));
		connect(HV,SIGNAL(activated()),this,SLOT(slot_histogram_view_activated()));
		d->m_histogram_views << HV;
	}
}

int MVCrossCorrelogramsWidget::currentUnit()
{
	return d->m_current_unit_num;
}

void MVCrossCorrelogramsWidget::setCurrentUnit(int num)
{
	if (d->m_current_unit_num==num) return;
	if (num<1) return;
	if (num>d->m_histogram_views.count()) return;

	d->m_current_unit_num=num;
	d->do_highlighting();
	emit currentUnitChanged();
}

int MVCrossCorrelogramsWidget::baseUnit()
{
	return d->m_base_unit_num;
}

void MVCrossCorrelogramsWidget::setBaseUnit(int num)
{
	d->m_base_unit_num=num;
}

void MVCrossCorrelogramsWidget::slot_histogram_view_clicked()
{
	int num=sender()->property("unit_number").toInt();
	setCurrentUnit(num);
}

void MVCrossCorrelogramsWidget::slot_histogram_view_activated()
{
	emit unitActivated(currentUnit());
}

void MVCrossCorrelogramsWidget::keyPressEvent(QKeyEvent *evt)
{
	if (evt->key()==Qt::Key_Left) {
		setCurrentUnit(this->currentUnit()-1);
	}
	else if (evt->key()==Qt::Key_Right) {
		setCurrentUnit(this->currentUnit()+1);
	}
	else if (evt->key()==Qt::Key_Up) {
		setCurrentUnit(this->currentUnit()-d->m_num_columns);
	}
	else if (evt->key()==Qt::Key_Down) {
		setCurrentUnit(this->currentUnit()+d->m_num_columns);
	}
}


void MVCrossCorrelogramsWidgetPrivate::do_highlighting()
{
	for (int i=0; i<m_histogram_views.count(); i++) {
		HistogramView *HV=m_histogram_views[i];
		int k=HV->property("unit_number").toInt();
		if (k==m_current_unit_num) {
			HV->setHighlighted(true);
		}
		else {
			HV->setHighlighted(false);
		}
	}
}
