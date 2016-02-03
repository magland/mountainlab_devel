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
#include "mvoverviewwidget.h"
#include "mvlabelcomparewidget.h"

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

    printf("testing\n");

	QStringList args;
	for (int i=1; i<argc; i++) {
		args << QString(argv[i]);
	}

	qsrand(QDateTime::currentDateTime().toMSecsSinceEpoch());

	QString working_path=CLP.named_parameters.value("working_path");
	QString output_path=CLP.named_parameters.value("output_path");

    QString mode=CLP.named_parameters.value("mode","overview");
    QString templates_path=CLP.named_parameters.value("templates");
    QString locations_path=CLP.named_parameters.value("locations");
    QString raw_path=CLP.named_parameters.value("raw");
    QString times_path=CLP.named_parameters.value("times");
    QString labels_path=CLP.named_parameters.value("labels");
    QString clusters_path=CLP.named_parameters.value("clusters");
    if (clusters_path.isEmpty()) clusters_path=CLP.named_parameters.value("cluster"); //historical compatibility
    QString primary_channels_path=CLP.named_parameters.value("primary-channels");
    QString cross_correlograms_path=CLP.named_parameters.value("cross-correlograms");
	if (cross_correlograms_path.isEmpty()) cross_correlograms_path=CLP.named_parameters.value("cross_correlograms");
	QString clips_path=CLP.named_parameters.value("clips");
	QString clips_index_path=CLP.named_parameters.value("clips-index");
	if (clips_index_path.isEmpty()) clips_index_path=CLP.named_parameters.value("clips_index");

    QString cluster2_path=CLP.named_parameters.value("cluster2"); //for mode=compare_labels

    if (mode=="overview") {
        MVOverviewWidget *W=new MVOverviewWidget;
        W->setWindowTitle(CLP.named_parameters.value("window_title","MountainView"));
        W->show();
        W->move(QApplication::desktop()->screen()->rect().topLeft()+QPoint(200,200));
        W->resize(1800,1200);
        if (!templates_path.isEmpty()) {
            Mda X; X.read(templates_path);
            W->setTemplates(X);
        }
        if (!locations_path.isEmpty()) {
            Mda X; X.read(locations_path);
            W->setElectrodeLocations(X);
        }
        if (!raw_path.isEmpty()) {
            DiskArrayModel *X=new DiskArrayModel;
            X->setPath(raw_path);
            W->setRaw(X,true);
        }
        if (!clips_path.isEmpty()) {
            DiskArrayModel *X=new DiskArrayModel;
            X->setPath(clips_path);
            W->setClips(X,true);
        }
        if (!clips_index_path.isEmpty()) {
            Mda X; X.read(clips_index_path);
            W->setClipsIndex(X);
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
            W->setTimesLabels(T,L);
        }
        if (!clusters_path.isEmpty()) {
            Mda CC; CC.read(clusters_path);
            int num_events=CC.N2();
            Mda T,L;
            T.allocate(1,num_events);
            L.allocate(1,num_events);
            for (int i=0; i<num_events; i++) {
                T.setValue(CC.value(1,i),0,i);
                L.setValue(CC.value(2,i),0,i);
            }
            W->setTimesLabels(T,L);
        }
        {
            Mda PC; PC.read(primary_channels_path);
            W->setPrimaryChannels(PC);
        }
        {
            W->setCrossCorrelogramsPath(cross_correlograms_path);
        }
        W->updateWidgets();
    }
    else if (mode=="compare_labels") {
        printf("compare_labels...\n");
        MVLabelCompareWidget *W=new MVLabelCompareWidget;
        W->setWindowTitle(CLP.named_parameters.value("window_title","MountainView - Compare Labels"));
        W->show();
        W->move(QApplication::desktop()->screen()->rect().topLeft()+QPoint(200,200));
        W->resize(1800,1200);
        if (!locations_path.isEmpty()) {
            Mda X; X.read(locations_path);
            W->setElectrodeLocations(X);
        }
        if (!raw_path.isEmpty()) {
            DiskArrayModel *X=new DiskArrayModel;
            X->setPath(raw_path);
            W->setRaw(X,true);
        }
        if ((!clusters_path.isEmpty())&&(!cluster2_path.isEmpty())) {
            Mda T1,L1,T2,L2;
            {
                Mda CC; CC.read(clusters_path);
                int num_events=CC.N2();
                Mda T,L;
                T.allocate(1,num_events);
                L.allocate(1,num_events);
                for (int i=0; i<num_events; i++) {
                    T.setValue(CC.value(1,i),0,i);
                    L.setValue(CC.value(2,i),0,i);
                }
                T1=T; L1=L;
            }
            {
                Mda CC; CC.read(cluster2_path);
                int num_events=CC.N2();
                Mda T,L;
                T.allocate(1,num_events);
                L.allocate(1,num_events);
                for (int i=0; i<num_events; i++) {
                    T.setValue(CC.value(1,i),0,i);
                    L.setValue(CC.value(2,i),0,i);
                }
                T2=T; L2=L;
            }
            W->setTimesLabels(T1,L1,T2,L2);
        }
        W->updateWidgets();
    }

	int ret=a.exec();

	printf("Number of files open: %d, number of unfreed mallocs: %d, number of unfreed megabytes: %g\n",jnumfilesopen(),jmalloccount(),(int)jbytesallocated()*1.0/1000000);

	return ret;
}
