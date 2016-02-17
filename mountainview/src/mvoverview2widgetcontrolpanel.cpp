#include "mvoverview2widgetcontrolpanel.h"

#include <QVBoxLayout>
#include <QGridLayout>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QMap>
#include <QCheckBox>

class MVOverview2WidgetControlPanelPrivate {
public:
	MVOverview2WidgetControlPanel *q;

	QMap<QString,QLineEdit *> m_lineedit_parameters;
    QMap<QString,QCheckBox *> m_checkbox_parameters;
    QMap<QString,QPushButton *> m_buttons;

	void add_group_label(QGridLayout *G,QString label);
    QCheckBox *add_check_box(QGridLayout *G,QString name,QString label,bool val);
	QLineEdit *add_int_box(QGridLayout *G,QString name,QString label,int val,int minval,int maxval);
	QLineEdit *add_float_box(QGridLayout *G,QString name,QString label,float val,float minval,float maxval);
	QPushButton *add_button(QGridLayout *G,QString name,QString label);
    void add_horizontal_divider(QVBoxLayout *layout);
};

void MVOverview2WidgetControlPanelPrivate::add_group_label(QGridLayout *G,QString label) {
	int r=G->rowCount();
	QLabel *X=new QLabel(label);
	QFont f=X->font();
	f.setPointSize(16);
	X->setFont(f);
    G->addWidget(X,r,0,1,2);
}

QCheckBox *MVOverview2WidgetControlPanelPrivate::add_check_box(QGridLayout *G, QString name, QString label, bool val)
{
    int r=G->rowCount();
    QCheckBox *X=new QCheckBox;
    X->setChecked(val);
    X->setText(QString("%1").arg(label));
    G->addWidget(X,r,1);
    m_checkbox_parameters[name]=X;
    X->setProperty("signal",name);
    q->connect(X,SIGNAL(toggled(bool)),q,SLOT(slot_checkbox_clicked(bool)));
    return X;
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
    m_buttons[name]=X;
    return X;
}

void MVOverview2WidgetControlPanelPrivate::add_horizontal_divider(QVBoxLayout *layout)
{
    QFrame *line = new QFrame;
    line->setFrameShape(QFrame::HLine); // Horizontal line
    line->setFrameShadow(QFrame::Sunken);
    line->setLineWidth(1);
    layout->addSpacing(25);
    layout->addWidget(line);
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
        d->add_horizontal_divider(layout);
	}


	{ // Templates
		QGridLayout *G=new QGridLayout;
        layout->addLayout(G);

		d->add_group_label(G,"Templates");
        d->add_int_box(G,"clip_size","Clip Size",100,20,10000)->setToolTip("Number of time points in clips");
		d->add_button(G,"update_templates","Update");
        d->add_horizontal_divider(layout);
	}

	{ // Amplitude Splitting
		QGridLayout *G=new QGridLayout;
        layout->addLayout(G);

		d->add_group_label(G,"Amplitude Splitting");
        d->add_check_box(G,"use_amplitude_split","Use amplitude split",false)->setToolTip("Split into peak amplitude shells.");
        d->add_float_box(G,"shell_width","Shell Width",1.5,0.1,20)->setToolTip("The width (in amplitude) of each shell");
        d->add_int_box(G,"min_per_shell","Min per shell",150,0,1500)->setToolTip("The minimum number of points in each shell");
        d->add_int_box(G,"min_amplitude","Min amplitude",0,0,100)->setToolTip("The minimum peak amplitude to include");
        d->add_button(G,"update_amplitude_split","Update");
        d->add_horizontal_divider(layout);
	}

	{ // Actions
        QGridLayout *G=new QGridLayout;
		layout->addLayout(G);

        d->add_group_label(G,"Actions");
        d->add_button(G,"open_templates","Open Templates")->setToolTip("Open a new window with auto-computed templates");
        d->add_button(G,"open_auto_correlograms","Open Auto-Correlograms")->setToolTip("Open a new auto-correlograms window");
        d->add_button(G,"open_raw_data","Open Raw Data")->setToolTip("Open a window of raw data");
        d->add_button(G,"open_clips","Open Clips")->setToolTip("Open clips for currently selected neuron.");
        d->add_horizontal_divider(layout);
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
    if (d->m_checkbox_parameters.contains(name)) return d->m_checkbox_parameters[name]->isChecked();
    return "";
}

void MVOverview2WidgetControlPanel::setParameterValue(QString name, QVariant val)
{
    if (d->m_lineedit_parameters.contains(name)) return d->m_lineedit_parameters[name]->setText(val.toString());
    if (d->m_checkbox_parameters.contains(name)) return d->m_checkbox_parameters[name]->setChecked(val.toBool());
}

void MVOverview2WidgetControlPanel::setParameterLabel(QString name, QString text)
{
    if (d->m_checkbox_parameters.contains(name)) return d->m_checkbox_parameters[name]->setText(text);
    if (d->m_buttons.contains(name)) return d->m_buttons[name]->setText(text);
}

void MVOverview2WidgetControlPanel::slot_button_clicked()
{
    emit signalButtonClicked(sender()->property("signal").toString());
}

void MVOverview2WidgetControlPanel::slot_checkbox_clicked(bool val)
{
    ((QCheckBox *)sender())->setChecked(val); //make sure this is set before we emit the signal
    emit signalButtonClicked(sender()->property("signal").toString());
}

