#include <QApplication>
#include <QScriptEngine>
#include <QDebug>
#include <qdatetime.h>
#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QProcess>
#include <QStringList>
#include "textfile.h"
#include "usagetracking.h"
#include "cvcommon.h"
#include "mountainviewwidget.h"
#include "mda.h"
#include <QDesktopWidget>
#include "get_command_line_params.h"
#include "diskarraymodel.h"
#include "histogramview.h"

/*
 * TO DO:
 * Clean up temporary files
 * */

bool download_file(QString url,QString fname) {
	QStringList args; args << "-o" << fname << url;
	int ret=QProcess::execute("/usr/bin/curl",args);
	if (ret!=0) {
		if (QFile::exists(fname)) {
			QFile::remove(fname);
		}
		return false;
	}
	return QFile::exists(fname);
}

void test_histogramview() {

	int N=100;
	float values[N];
	for (int i=0; i<N; i++) {
		values[i]=(qrand()%10000)*1.0/10000;
	}

	HistogramView *W=new HistogramView;
	W->setData(N,values);
	W->autoSetBins(N/5);
	W->show();
}


int main(int argc, char *argv[]) {
	QApplication a(argc, argv);
	//MainWindow w;
	//w.show();

//	{
//		test_histogramview();
//		return a.exec();
//	}


    QStringList required;
    QStringList optional;
    CLParams CLP=get_command_line_params(argc,argv,required);

	QStringList args;
	for (int i=1; i<argc; i++) {
		args << QString(argv[i]);
	}

	qsrand(QDateTime::currentDateTime().toMSecsSinceEpoch());

	QString working_path=CLP.named_parameters.value("working_path");
	QString output_path=CLP.named_parameters.value("output_path");

    qDebug() << CLP.named_parameters;

    QString templates_path=CLP.named_parameters.value("templates");
    QString templates_whitened_path=CLP.named_parameters.value("templates-whitened");
    QString locations_path=CLP.named_parameters.value("locations");
    QString raw_path=CLP.named_parameters.value("raw");
    QString raw_whitened_path=CLP.named_parameters.value("raw-whitened");
    QString times_path=CLP.named_parameters.value("times");
    QString labels_path=CLP.named_parameters.value("labels");
    QString cluster_path=CLP.named_parameters.value("cluster");
    QString primary_channels_path=CLP.named_parameters.value("primary-channels");
    QString cross_correlograms_path=CLP.named_parameters.value("cross-correlograms");

    if (!QFile::exists(templates_whitened_path)) templates_whitened_path="";
    if (!QFile::exists(raw_whitened_path)) raw_whitened_path="";

    qDebug() << templates_path;

	MountainViewWidget W;
    W.show();
    W.move(QApplication::desktop()->screen()->rect().topLeft()+QPoint(300,300));
    if (!templates_path.isEmpty()) {
        Mda X; X.read(templates_path);
        W.setTemplates(X);
    }
    if (!templates_whitened_path.isEmpty()) {
        Mda X; X.read(templates_whitened_path);
        W.setTemplatesWhitened(X);
    }
    if (!locations_path.isEmpty()) {
        Mda X; X.read(locations_path);
		W.setElectrodeLocations(X);
    }
    if (!raw_path.isEmpty()) {
        DiskArrayModel *X=new DiskArrayModel;
        X->setPath(raw_path);
        X->createFileHierarchyIfNeeded();
        W.setRaw(X);
    }
    if (!raw_whitened_path.isEmpty()) {
        DiskArrayModel *X=new DiskArrayModel;
        X->setPath(raw_whitened_path);
        X->createFileHierarchyIfNeeded();
        W.setRawWhitened(X);
    }
    if (!times_path.isEmpty()) {
        Mda T; T.read(times_path);
        Mda L;
        if (!labels_path.isEmpty()) {
            L.read(labels_path);
        }
        else {
            L.allocate(T.N1(),T.N2());
            for (int ii=0; ii<L.totalSize(); ii++) L.setValue1(1,ii);
        }
        W.setTimesLabels(T,L);
    }
    if (!cluster_path.isEmpty()) {
        Mda CC; CC.read(cluster_path);
        int num_events=CC.N2();
        Mda T,L;
        T.allocate(1,num_events);
        L.allocate(1,num_events);
        for (int i=0; i<num_events; i++) {
            T.setValue(CC.value(1,i),0,i);
            L.setValue(CC.value(2,i),0,i);
        }
        W.setTimesLabels(T,L);
    }
	{
		Mda PC; PC.read(primary_channels_path);
		W.setPrimaryChannels(PC);
	}
    {
        W.setCrossCorrelogramsPath(cross_correlograms_path);
    }


	int ret=a.exec();

	printf("Number of files open: %d, number of unfreed mallocs: %d, number of unfreed megabytes: %g\n",jnumfilesopen(),jmalloccount(),(int)jbytesallocated()*1.0/1000000);

	return ret;
}
