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


int main(int argc, char *argv[]) {
	QApplication a(argc, argv);
	//MainWindow w;
	//w.show();

    QStringList required;
	QStringList optional; optional << "working_path" << "output_path";
    CLParams CLP=get_command_line_params(argc,argv,required,optional);

	QStringList args;
	for (int i=1; i<argc; i++) {
		args << QString(argv[i]);
	}

	qsrand(QDateTime::currentDateTime().toMSecsSinceEpoch());

	QString working_path=CLP.named_parameters.value("working_path","/home/magland/dev/ms11d45A/working");
	QString output_path=CLP.named_parameters.value("output_path","/home/magland/dev/ms11d45A/output");

	QString templates_path=QString("%1/templates.mda").arg(output_path);
	QString locations_path=QString("%1/locations.mda").arg(output_path);
	QString raw_path=QString("%1/raw.mda").arg(output_path);
	QString times_path=QString("%1/times.mda").arg(output_path);
	QString labels_path=QString("%1/labels.mda").arg(output_path);
	QString primary_channels_path=QString("%1/primary_channels.mda").arg(output_path);

	MountainViewWidget W;
    W.show();
    W.move(QApplication::desktop()->screen()->rect().topLeft()+QPoint(300,300));
    if (!templates_path.isEmpty()) {
        Mda X; X.read(templates_path);
        W.setTemplates(X);
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
	{
		Mda PC; PC.read(primary_channels_path);
		W.setPrimaryChannels(PC);
	}


	int ret=a.exec();

	printf("Number of files open: %d, number of unfreed mallocs: %d, number of unfreed megabytes: %g\n",jnumfilesopen(),jmalloccount(),(int)jbytesallocated()*1.0/1000000);

	return ret;
}
