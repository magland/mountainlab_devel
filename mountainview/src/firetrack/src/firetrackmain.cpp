#include <QApplication>
#include <QScriptEngine>
#include <QDebug>
#include <qdatetime.h>
#include "ftcontroller.h"
#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QProcess>
#include <QStringList>
#include "textfile.h"
#include "usagetracking.h"
#include "cvcommon.h"
#include "ftelectrodearrayview.h"

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

	QStringList args;
	for (int i=1; i<argc; i++) {
		args << QString(argv[i]);
	}

	qsrand(QDateTime::currentDateTime().toMSecsSinceEpoch());

	//args << "/tmp/tp903426a1_f13b_43a2_a1b8_5ff229f26018_run_java_script.js";

//	FTElectrodeArrayView *V=new FTElectrodeArrayView;
//	int M=100,T=30;
//	Mda L; L.allocate(M,2);
//	Mda W; W.allocate(M,T);
//	for (int i=0; i<M; i++) {
//		L.setValue(i/6 + ((i%6)%2)*0.5,i,0);
//		L.setValue(i%6,i,1);
//		for (int t=0; t<T; t++) {
//			W.setValue(i+t,i,t);
//		}
//	}
//	V->setElectrodeLocations(L);
//	V->setWaveform(W);
//	V->show();
//	return a.exec();

	QString script_path;
	QString waveforms_path;
	QString locations_path;
	for (int i=0; i<args.count(); i++) {
		QString str=args[i];
		if (str=="--waveforms") {
			waveforms_path=args.value(i+1);
		}
		if (str=="--locations") {
			locations_path=args.value(i+1);
		}
	}

	QScriptEngine *engine=new QScriptEngine;

	//qScriptRegisterMetaType(engine, myObjectToScriptValue, myObjectFromScriptValue);


	FTController FIRETRACK;
	QScriptValue FIRETRACK_value = engine->newQObject(&FIRETRACK);
	   //QScriptValue FIRETRACK_value = engine->newQObject(new QObject());
	engine->globalObject().setProperty("FIRETRACK", FIRETRACK_value);

	QString script;
	if (!script_path.isEmpty()) {
		script=read_text_file(script_path);
	} else {
		if (waveforms_path.isEmpty()) {
            waveforms_path=a.applicationDirPath()+"/../testdata/waveforms_first_5e5_points.mda";
            if ((!QFile::exists(waveforms_path))||(QFileInfo(waveforms_path).size()<1e6)) {
				QMessageBox::StandardButton reply=
				QMessageBox::question(0,"Start download?","Test data must be downloaded. Start download?",QMessageBox::Yes|QMessageBox::No,QMessageBox::Yes);
				if (reply==QMessageBox::No) exit(0);
				QMessageBox::information(0,"Downloading","Please wait while file is downloaded. Click OK to start download.");
				QString url="http://97.107.129.125/waveforms_first_5e5_points.mda";
				if (!download_file(url,waveforms_path)) {
					QMessageBox::information(0,"Error Downloading","There was a problem downloading the file "+url+". Is curl installed? Exiting.");
					exit(0);
				}
				QMessageBox::information(0,"File downloaded","It appears that the file was downloaded properly. Press OK to continue.");
			}

            //waveforms_path="/home/magland/gazelle_backup/current/matlab/scda_ss/jfm/core/scratch/test_all_channels/first_1e3_points_filtered.mda";
		}
		if (locations_path.isEmpty()) {
			locations_path=a.applicationDirPath()+"/../testdata/locations.mda";
		}
		script="";
		script+=QString("var V%1=FIRETRACK.createFireTrackWidget();\n").arg(0);
		script+=QString("var X%1=FIRETRACK.readArray('%2');\n").arg(0).arg(waveforms_path);
		script+=QString("var L%1=FIRETRACK.readArray('%2');\n").arg(0).arg(locations_path);
		script+=QString("V%1.setWaveforms(X%1);\n").arg(0);
		script+=QString("V%1.setElectrodeLocations(L%1);\n").arg(0);
		//script+=QString("V%1.setTitle(\"This is a test %1.\");\n").arg(0);
		script+=QString("V%1.show();\n").arg(0);
        script+=QString("V%1.animate();\n").arg(0);
		script+="\n";
	}

	QScriptValue result = engine->evaluate(script);
	if (result.isError()) {
		qWarning() << "Error running script: "+result.toString();
	}

	CleanupObject cleanup_object;
	QObject::connect(&a, SIGNAL(aboutToQuit()), &cleanup_object, SLOT(closing()));

	int ret=a.exec();

	engine->collectGarbage();

	delete engine;
	printf("Number of files open: %d, number of unfreed mallocs: %d, number of unfreed megabytes: %g\n",jnumfilesopen(),jmalloccount(),(int)jbytesallocated()*1.0/1000000);

	return ret;
}
