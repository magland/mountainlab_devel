#include <stdlib.h>
#include <stdio.h>

#include "get_command_line_params.h"
#include "processtracker.h"
#include <QCoreApplication>

#include "extract.h"
#include "bandpass_filter.h"
#include "normalize_channels.h"
#include "whiten.h"
#include "detect.h"
#include "features0.h"
#include "cluster.h"

void register_processors(ProcessTracker &PT) {
	{
		PTProcessor P;
		P.command="extract";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.1";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="bandpass_filter";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.1";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="normalize_channels";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.1";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="whiten";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.1";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="detect";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.11";
		PT.registerProcessor(P);
	}
    {
        PTProcessor P;
        P.command="features";
        P.input_file_pnames << "input";
        P.input_file_pnames << "detect";
        P.output_file_pnames << "output";
        P.version="0.1";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="cluster";
        P.input_file_pnames << "input";
        P.output_file_pnames << "output";
        P.version="0.1";
        PT.registerProcessor(P);
    }
}

void extract_usage() {
	printf("mountainsort extract --input=in.mda --output=out.mda --num_channels=72 --channels=1,2,3,4 --t1=0 --t2=1e6\n");
}

void bandpass_filter_usage() {
	printf("mountainsort bandpass_filter --input=in.mda --output=out.mda --samplefreq=30000 --freq_min=300 --freq_max=10000 --outlier_threshold=0\n");
}

void normalize_channels_usage() {
	printf("mountainsort normalize_channels --input=in.mda --output=out.mda\n");
}

void whiten_usage() {
	printf("mountainsort whiten --input=in.mda --output=out.mda --ncomp=4\n");
}

void detect_usage() {
	printf("mountainsort detect --input=in.mda --output=out.mda --inner_window_width=40 --outer_window_width=1000 --threshold=5\n");
}

void features_usage() {
    printf("mountainsort features --input=in.mda --detect=detect.mda --output=out.mda --num_features=6 --clip_size=100\n");
}

void cluster_usage() {
    printf("mountainsort cluster --input=in.mda --output=out.mda\n");
}

int main(int argc,char *argv[]) {

	QCoreApplication app(argc,argv); //important for qApp->applicationDirPath() in processtracker

	CLParams CLP;
	QStringList required;
	CLP=get_command_line_params(argc,argv,required);

	ProcessTracker PT;
	register_processors(PT);

	if (CLP.unnamed_parameters.count()>1) {
		printf("Only one command parameter may be specified.\n");
		return -1;
	}

	QString command=CLP.unnamed_parameters.value(0);

	if (command.isEmpty()) {
		printf("\nmountainsort processors:\n");
		for (int i=0; i<PT.processorCount(); i++) {
			QString cmd=PT.processor(i).command;
			printf("%s ",cmd.toLatin1().data());
		}
		printf("\n\n");
		return -1;
	}

	PTProcessor PP=PT.findProcessor(command);
	printf("%s version %s\n",PP.command.toLatin1().data(),PP.version.toLatin1().data());

	if (!CLP.named_parameters.value("force").toInt()) {
		if (PT.processAlreadyCompleted(CLP)) {
			printf("Process already completed.\n");
			return 0;
		}
	}
	CLP.named_parameters.remove("force");

	if (command=="extract") {
		QString input_path=CLP.named_parameters["input"];
		QString output_path=CLP.named_parameters["output"];
		int num_channels=CLP.named_parameters["num_channels"].toInt();
		long t1=CLP.named_parameters["t1"].toLong();
		long t2=CLP.named_parameters["t2"].toLong();
		QStringList channels_str=CLP.named_parameters["channels"].split(",");
		int M=channels_str.count();
		int channels[M];
		for (int m=0; m<M; m++) channels[m]=channels_str[m].toInt();

		if ((input_path.isEmpty())||(output_path.isEmpty())) {extract_usage(); return -1;}
		if (M==0) {extract_usage(); return -1;}

		if (!extract(input_path.toLatin1().data(),output_path.toLatin1().data(),num_channels,M,channels,t1,t2)) {
			printf("Error in extract.\n");
			return -1;
		}
	}
	else if (command=="bandpass_filter") {
		QString input_path=CLP.named_parameters["input"];
		QString output_path=CLP.named_parameters["output"];
		double samplefreq=CLP.named_parameters["samplefreq"].toDouble();
		double freq_min=CLP.named_parameters["freq_min"].toDouble();
		double freq_max=CLP.named_parameters["freq_max"].toDouble();
		double outlier_threshold=CLP.named_parameters["outlier_threshold"].toDouble();

		if ((input_path.isEmpty())||(output_path.isEmpty())) {bandpass_filter_usage(); return -1;}
		if ((samplefreq==0)||(freq_min==0)||(freq_max==0)) {bandpass_filter_usage(); return -1;}

		if (!bandpass_filter(input_path.toLatin1().data(),output_path.toLatin1().data(),samplefreq,freq_min,freq_max,outlier_threshold)) {
			printf("Error in bandpass_filter.\n");
			return -1;
		}
	}
	else if (command=="normalize_channels") {
		QString input_path=CLP.named_parameters["input"];
		QString output_path=CLP.named_parameters["output"];

		if ((input_path.isEmpty())||(output_path.isEmpty())) {normalize_channels_usage(); return -1;}

		if (!normalize_channels(input_path.toLatin1().data(),output_path.toLatin1().data())) {
			printf("Error in normalize_channels.\n");
			return -1;
		}
	}
	else if (command=="whiten") {
		QString input_path=CLP.named_parameters["input"];
		QString output_path=CLP.named_parameters["output"];
		int ncomp=CLP.named_parameters["ncomp"].toInt();

		if ((input_path.isEmpty())||(output_path.isEmpty())) {whiten_usage(); return -1;}
		if (ncomp==0) {whiten_usage(); return -1;}

		if (!whiten(input_path.toLatin1().data(),output_path.toLatin1().data(),ncomp)) {
			printf("Error in whiten.\n");
			return -1;
		}
	}
	else if (command=="detect") {
		QString input_path=CLP.named_parameters["input"];
		QString output_path=CLP.named_parameters["output"];
		int inner_window_width=CLP.named_parameters["inner_window_width"].toInt();
		int outer_window_width=CLP.named_parameters["outer_window_width"].toInt();
		float threshold=CLP.named_parameters["threshold"].toFloat();

		if ((input_path.isEmpty())||(output_path.isEmpty())) {detect_usage(); return -1;}
		if (inner_window_width==0) {detect_usage(); return -1;}
		if (outer_window_width==0) {detect_usage(); return -1;}
		if (threshold==0) {detect_usage(); return -1;}

		if (!detect(input_path.toLatin1().data(),output_path.toLatin1().data(),inner_window_width,outer_window_width,threshold)) {
			printf("Error in detect.\n");
			return -1;
		}
	}
    else if (command=="features") {
        QString input_path=CLP.named_parameters["input"];
        QString detect_path=CLP.named_parameters["detect"];
        QString output_path=CLP.named_parameters["output"];
        int num_features=CLP.named_parameters["num_features"].toInt();
        int clip_size=CLP.named_parameters["clip_size"].toInt();

        if ((input_path.isEmpty())||(detect_path.isEmpty())||(output_path.isEmpty())) {features_usage(); return -1;}
        if (num_features==0) {features_usage(); return -1;}
        if (clip_size==0) {features_usage(); return -1;}

        if (!features(input_path.toLatin1().data(),detect_path.toLatin1().data(),output_path.toLatin1().data(),num_features,clip_size)) {
            printf("Error in features.\n");
            return -1;
        }
    }
    else if (command=="cluster") {
        QString input_path=CLP.named_parameters["input"];
        QString output_path=CLP.named_parameters["output"];

        if ((input_path.isEmpty())||(output_path.isEmpty())) {cluster_usage(); return -1;}

        if (!cluster(input_path.toLatin1().data(),output_path.toLatin1().data())) {
            printf("Error in cluster.\n");
            return -1;
        }
    }
	else {
		printf("Unknown command: %s\n",command.toLatin1().data());
		return -1;
	}

	PT.reportProcessCompleted(CLP);

	return 0;
}
