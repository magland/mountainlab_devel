#include "mvoverview2widgetcontrolpanel.h"

#include <QVBoxLayout>
#include <QGridLayout>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QMap>

class MVOverview2WidgetControlPanelPrivate {
public:
	MVOverview2WidgetControlPanel *q;

	QMap<QString,QLineEdit *> m_lineedit_parameters;

	void add_group_label(QGridLayout *G,QString label);
	QLineEdit *add_int_box(QGridLayout *G,QString name,QString label,int val,int minval,int maxval);
	QLineEdit *add_float_box(QGridLayout *G,QString name,QString label,float val,float minval,float maxval);
	QPushButton *add_button(QGridLayout *G,QString name,QString label);
};

void MVOverview2WidgetControlPanelPrivate::add_group_label(QGridLayout *G,QString label) {
	int r=G->rowCount();
	QLabel *X=new QLabel(label);
	QFont f=X->font();
	f.setPointSize(16);
	X->setFont(f);
	G->addWidget(X,r,0,1,2);
}

QLineEdit *MVOverview2WidgetControlPanelPrivate::add_int_box(QGridLayout *G,QString name,QString label,int val,int minval,int maxval) {
	Q_UNUSED(minval)
	Q_UNUSED(maxval)
	int r=G->rowCount();
	QLineEdit *X=new QLineEdit;
	X->setText(QString("%1").arg(val));
	G->addWidget(new QLabel(label),r,0);
	G->addWidget(X,r,1);
	m_lineedit_parameters[name]=X;
	return X;
}

QLineEdit *MVOverview2WidgetControlPanelPrivate::add_float_box(QGridLayout *G, QString name, QString label, float val, float minval, float maxval)
{
	Q_UNUSED(minval)
	Q_UNUSED(maxval)
	int r=G->rowCount();
	QLineEdit *X=new QLineEdit;
	X->setText(QString("%1").arg(val));
	G->addWidget(new QLabel(label),r,0);
	G->addWidget(X,r,1);
	m_lineedit_parameters[name]=X;
	return X;
}

QPushButton *MVOverview2WidgetControlPanelPrivate::add_button(QGridLayout *G,QString name,QString label) {
	int r=G->rowCount();
	QPushButton *X=new QPushButton(label);
	G->addWidget(X,r,1);
	X->setProperty("signal",name);
	q->connect(X,SIGNAL(clicked(bool)),q,SLOT(slot_button_clicked()));
	return X;
}

MVOverview2WidgetControlPanel::MVOverview2WidgetControlPanel(QWidget *parent) : QWidget(parent)
{
	d=new MVOverview2WidgetControlPanelPrivate;
	d->q=this;

	QVBoxLayout *layout=new QVBoxLayout;

	{ // Cross-correlograms
		QGridLayout *G=new QGridLayout;
		layout->addLayout(G);

		d->add_group_label(G,"Cross-correlograms");
		d->add_int_box(G,"max_dt","Max. dt (ms)",100,1,3000);
		d->add_button(G,"update_cross_correlograms","Update");
	}

	{ // Templates
		QGridLayout *G=new QGridLayout;
		layout->addLayout(G);

		d->add_group_label(G,"Templates");
		d->add_int_box(G,"clip_size","Clip Size",100,20,10000);
		d->add_button(G,"update_templates","Update");
	}

	{ // Amplitude Splitting
		QGridLayout *G=new QGridLayout;
		layout->addLayout(G);

		d->add_group_label(G,"Amplitude Splitting");
		d->add_float_box(G,"shell_width","Shell Width",1.5,0.1,20);
		d->add_int_box(G,"min_per_shell","Min per shell",150,0,1500);
		d->add_button(G,"amplitude_split","Split");
	}

	layout->addStretch(0);

	this->setLayout(layout);
}

MVOverview2WidgetControlPanel::~MVOverview2WidgetControlPanel()
{
	delete d;
}

QVariant MVOverview2WidgetControlPanel::getParameterValue(QString name)
{
	if (d->m_lineedit_parameters.contains(name)) return d->m_lineedit_parameters[name]->text();
	return "";
}

void MVOverview2WidgetControlPanel::slot_button_clicked()
{
	emit signalButtonClicked(sender()->property("signal").toString());
}

