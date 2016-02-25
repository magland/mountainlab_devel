#include <stdlib.h>
#include <stdio.h>
#include <QFile>

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
#include "isobranch.h"
#include "split_firings.h"
#include "templates.h"
#include "consolidate.h"
#include "extract_clips.h"
#include "extract_channels.h"
#include "remove_artifacts.h"
#include "create_clips_file.h"
#include "assemble_firings_file.h"
#include "get_principal_components.h"
#include "fit.h"
#include "cross_correlograms.h"
#include "confusion_matrix.h"
#include <QTime>
#include "Eigen/Core"
#include "Eigen/SVD"
#include "process_msh.h"
#include "isosplit2.h"
#include "branch_cluster_v1.h"

void register_processors(ProcessTracker &PT) {
	{
		PTProcessor P;
		P.command="extract";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.11";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="bandpass_filter";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.11";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="normalize_channels";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
		P.version="0.11";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="whiten";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
        P.version="0.16";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="detect";
		P.input_file_pnames << "input";
		P.output_file_pnames << "output";
        P.version="0.18";
		PT.registerProcessor(P);
	}
    {
        PTProcessor P;
        P.command="features";
        P.input_file_pnames << "input";
        P.input_file_pnames << "detect";
        P.input_file_pnames << "adjacency";
        P.output_file_pnames << "output";
		P.version="0.20";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="cluster";
        P.input_file_pnames << "input";
        P.output_file_pnames << "output";
		P.version="0.17";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="isobranch";
        P.input_file_pnames << "input_clips";
        P.output_file_pnames << "output_labels";
        P.version="0.12";
        PT.registerProcessor(P);
    }
	{
		PTProcessor P;
		P.command="split_firings";
		P.input_file_pnames << "input";
		P.input_file_pnames << "firings";
		P.output_file_pnames << "output";
		P.version="0.15";
		PT.registerProcessor(P);
	}
    {
        PTProcessor P;
        P.command="templates";
        P.input_file_pnames << "input";
		P.input_file_pnames << "firings";
        P.output_file_pnames << "output";
		P.version="0.17";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="consolidate";
		P.input_file_pnames << "firings";
        P.input_file_pnames << "templates";
        P.output_file_pnames << "cluster_out";
        P.output_file_pnames << "templates_out";
        P.output_file_pnames << "load_channels_out";
		P.version="0.47";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="fit";
        P.input_file_pnames << "input";
		P.input_file_pnames << "firings";
        P.output_file_pnames << "templates";
        P.output_file_pnames << "cluster_out";
		P.version="0.14";
        PT.registerProcessor(P);
    }
	{
		PTProcessor P;
		P.command="extract_clips";
		P.input_file_pnames << "input";
        P.input_file_pnames << "detect";
        P.output_file_pnames << "output_clips";
        P.version="0.15";
		PT.registerProcessor(P);
	}
    {
        PTProcessor P;
        P.command="extract_channels";
        P.input_file_pnames << "input";
        P.output_file_pnames << "output";
        P.version="0.10";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="remove_artifacts";
        P.input_file_pnames << "input";
        P.output_file_pnames << "output";
        P.version="0.10";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="create_clips_file";
        P.input_file_pnames << "input";
		P.input_file_pnames << "firings";
        P.output_file_pnames << "output";
        P.output_file_pnames << "index_out";
        P.version="0.1";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
		P.command="assemble_firings_file";
        P.input_file_pnames << "input_detect";
        P.input_file_pnames << "input_labels";
		P.output_file_pnames << "output_firings";
        P.version="0.11";
        PT.registerProcessor(P);
    }
	{
		PTProcessor P;
		P.command="cross_correlograms";
		P.input_file_pnames << "firings";
		P.output_file_pnames << "output";
		P.version="0.11";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="confusion_matrix";
		P.input_file_pnames << "firings1" << "firings2";
		P.output_file_pnames << "output";
		P.version="0.12";
		PT.registerProcessor(P);
	}
    {
        PTProcessor P;
        P.command="branch_cluster_v1";
        P.input_file_pnames << "raw" << "detect";
        P.output_file_pnames << "firings";
        P.version="0.1";
        PT.registerProcessor(P);
    }
	{
		PTProcessor P;
		P.command="copy";
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
	//printf("mountainsort whiten --input=in.mda --output=out.mda --ncomp=4\n");
	printf("mountainsort whiten --input=in.mda --output=out.mda\n");
}

void detect_usage() {
    printf("mountainsort detect --input=in.mda --output=out.mda --inner_window_width=40 --outer_window_width=1000 --threshold=3 --normalize=0 --individual_channels=1\n");
}

void features_usage() {
    printf("mountainsort features --input=in.mda --detect=detect.mda --adjacency=adjacency.mda --output=out.mda --num_features=6 --clip_size=100\n");
}

void cluster_usage() {
	printf("mountainsort cluster --input=in.mda --output=out.mda --ks_threshold=1.4 --K_init=25 \n");
}

void isobranch_usage() {
    printf("mountainsort isobranch --input_clips=clips.mda --output_labels=labels.mda --branch_thresholds=2.5,3,3.5,4,5 --isocut_threshold=1.2 --K_init=30 --num_features=3 \n");
}

void split_firings_usage() {
	printf("mountainsort split_firings --input=in.mda --cluster=cluster.mda --output=out.mda --num_features=3 --clip_size=100 --ks_threshold=1.4 --K_init \n");
}

void templates_usage() {
    printf("mountainsort templates --input=in.mda --cluster=cluster.mda --output=out.mda --clip_size=100\n");
}

void consolidate_usage() {
	printf("mountainsort consolidate --cluster=cluster.mda --templates=templates.mda --cluster_out=cluster_out.mda --templates_out=templates_out.mda --load_channels_out=load_channels.mda --coincidence_threshold=0.5\n");
}

void fit_usage() {
    printf("mountainsort fit --input=input.mda --cluster=cluster.mda --templates=templates.mda --cluster_out=cluster_out.mda\n");
}

void extract_clips_usage() {
    printf("mountainsort extract_clips --input=raw.mda --detect=detect.mda --output=clips.mda --clip_size=100\n");
}

void extract_channels_usage() {
    printf("mountainsort extract_channels --input=raw.mda --output=raw2.mda --channels=1,2,3,4,5\n");
}

void remove_artifacts_usage() {
    printf("mountainsort remove_artifacts --input=pre1a.mda --output=pre1.mda --threshold=8 --normalize=1 --exclude_interval=1000\n");
}

void create_clips_file_usage() {
    printf("mountainsort create_clips_file --input=raw.mda --cluster=cluster.mda --output=clips.mda --index_out=clips_index.mda --clip_size=100\n");
}

void assemble_firings_file_usage() {
	printf("mountainsort assemble_firings_file --input_detect=detect.mda --input_labels=labels.mda --output_firings=firings.mda\n");
}

void cross_correlograms_usage() {
	printf("mountainsort cross_correlograms --firings=firings.mda --output=cross_correlograms.mda --max_dt=1500\n");
}

void confusion_matrix_usage() {
	printf("mountainsort confusion_matrix --firings1=firings1.mda --firings2=firings2.mda --output=confusion_matrix.mda --max_matching_offset=3\n");
}

void branch_cluster_v1_usage() {
    printf("mountainsort branch_cluster_v1 --raw=pre.mda --detect=detect.mda --adjacency_matrix= --firings=firings.mda --clip_size=100 --min_section_count=50 --section_increment=0.5 --num_features=6\n");
}

void copy_usage() {
	printf("mountainsort copy --input=in.mda --output=out.mda \n");
}

void test_svd() {
	/*
	 Compare with this MATLAB script:
		m=zeros(3,5);
		for r=1:size(m,1)
			for c=1:size(m,2)
				if (r==c) m(r,c)=5;
				else m(r,c)=r+c;
				end;
			end;
		end;

		[U,D,V]=svd(m,'econ');
		diag(D)
		U
		V
	*/


	using namespace Eigen;
	MatrixXf m(3,5);
	int ii=0;
	for (int c=0; c<m.cols(); c++) {
		for (int r=0; r<m.rows(); r++) {
			if (r==c) {
				m.data()[ii]=5;
			}
			else {
				m.data()[ii]=(r+1)+(c+1);
			}
			ii++;
		}
	}
	//cout << "Here is the matrix m:" << endl << m << endl;
	JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
	MatrixXf sv=svd.singularValues();
	printf("The singular values are %ld x %ld --- %g,%g,%g\n",sv.rows(),sv.cols(),sv.data()[0],sv.data()[1],sv.data()[2]);
	MatrixXf U=svd.matrixU();
	printf("The matrix U is %ld x %ld\n",U.rows(),U.cols());
	MatrixXf V=svd.matrixV();
	printf("The matrix V is %ld x %ld\n",V.rows(),V.cols());


	printf("\n\nU:\n");
	for (int r=0; r<U.rows(); r++) {
		for (int c=0; c<U.cols(); c++) {
			printf("%g  ",U.data()[r+U.rows()*c]);
		}
		printf("\n");
	}

	printf("\n\nV:\n");
	for (int r=0; r<V.rows(); r++) {
		for (int c=0; c<V.cols(); c++) {
			printf("%g  ",V.data()[r+V.rows()*c]);
		}
		printf("\n");
	}
}

int main(int argc,char *argv[]) {

	QCoreApplication app(argc,argv); //important for qApp->applicationDirPath() in processtracker

    //int test_max_int=1e9;
    //printf("%ld,%ld\n",test_max_int+1,sizeof(test_max_int));
    //return 0;

	//test_svd();

	CLParams CLP;
	QStringList required;
	CLP=get_command_line_params(argc,argv,required);

    if (CLP.unnamed_parameters.value(0)=="test_isosplit2_routines") {
        test_isosplit2_routines();
        return 0;
    }

	ProcessTracker PT;
	register_processors(PT);

    if (CLP.unnamed_parameters.value(0).endsWith(".msh")) {
        return process_msh(CLP.unnamed_parameters.value(0),argc,argv);
    }

	if (CLP.unnamed_parameters.count()>1) {
		printf("Only one command parameter may be specified.\n");
		qDebug()  << CLP.unnamed_parameters;
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

	QTime timer; timer.start();

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
		//int ncomp=CLP.named_parameters["ncomp"].toInt();

		if ((input_path.isEmpty())||(output_path.isEmpty())) {whiten_usage(); return -1;}
		//if (ncomp==0) {whiten_usage(); return -1;}

		//if (!whiten(input_path.toLatin1().data(),output_path.toLatin1().data(),ncomp)) {
		if (!whiten(input_path.toLatin1().data(),output_path.toLatin1().data())) {
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
        int normalize=CLP.named_parameters["normalize"].toInt();
		int individual_channels=1;
		if (CLP.named_parameters.contains("individual_channels"))
			individual_channels=CLP.named_parameters["individual_channels"].toInt();

		if ((input_path.isEmpty())||(output_path.isEmpty())) {detect_usage(); return -1;}
		if (inner_window_width==0) {detect_usage(); return -1;}
		if (outer_window_width==0) {detect_usage(); return -1;}
		if (threshold==0) {detect_usage(); return -1;}

        if (!detect(input_path.toLatin1().data(),output_path.toLatin1().data(),inner_window_width,outer_window_width,threshold,normalize,(individual_channels!=0))) {
			printf("Error in detect.\n");
			return -1;
		}
	}
    else if (command=="features") {
        QString input_path=CLP.named_parameters["input"];
        QString detect_path=CLP.named_parameters["detect"];
        QString adjacency_path=CLP.named_parameters["adjacency"];
        QString output_path=CLP.named_parameters["output"];
        int num_features=CLP.named_parameters["num_features"].toInt();
        int clip_size=CLP.named_parameters["clip_size"].toInt();

        if ((input_path.isEmpty())||(detect_path.isEmpty())||(adjacency_path.isEmpty())||(output_path.isEmpty())) {features_usage(); return -1;}
        if (num_features==0) {features_usage(); return -1;}
        if (clip_size==0) {features_usage(); return -1;}

        if (!features(input_path.toLatin1().data(),detect_path.toLatin1().data(),adjacency_path.toLatin1().data(),output_path.toLatin1().data(),num_features,clip_size)) {
            printf("Error in features.\n");
            return -1;
        }
    }
    else if (command=="cluster") {
        QString input_path=CLP.named_parameters["input"];
        QString output_path=CLP.named_parameters["output"];

        if ((input_path.isEmpty())||(output_path.isEmpty())) {cluster_usage(); return -1;}
		float ks_threshold=CLP.named_parameters["ks_threshold"].toFloat();
		int K_init=CLP.named_parameters["K_init"].toInt();
		if (ks_threshold==0) ks_threshold=1.4;
		if (K_init==0) K_init=25;

		if (!cluster(input_path.toLatin1().data(),output_path.toLatin1().data(),ks_threshold,K_init)) {
            printf("Error in cluster.\n");
            return -1;
        }
    }
    else if (command=="isobranch") {
        QString input_clips_path=CLP.named_parameters["input_clips"];
        QString output_labels_path=CLP.named_parameters["output_labels"];

        if ((input_clips_path.isEmpty())||(output_labels_path.isEmpty())) {isobranch_usage(); return -1;}
        float isocut_threshold=CLP.named_parameters["isocut_threshold"].toFloat();
        int K_init=CLP.named_parameters["K_init"].toInt();
        int num_features=CLP.named_parameters["num_features"].toInt();
        QList<QString> branch_thresholds_str=CLP.named_parameters["branch_thresholds"].split(",");
        QList<float> branch_thresholds;
        for (int i=0; i<branch_thresholds_str.count(); i++) {
            branch_thresholds << (branch_thresholds_str[i].trimmed()).toFloat();
        }
        if ((isocut_threshold==0)||(K_init==0)||(num_features==0)||(branch_thresholds.count()==0)) {
            isobranch_usage(); return -1;
        }

        if (!isobranch(input_clips_path.toLatin1(),output_labels_path.toLatin1(),branch_thresholds,num_features,isocut_threshold,K_init)) {
            printf("Error in isobranch.\n");
            return -1;
        }
    }
	else if (command=="split_firings") {
		QString input_path=CLP.named_parameters["input"];
		QString firings_path=CLP.named_parameters["firings"];
		QString output_path=CLP.named_parameters["output"];
		int num_features=CLP.named_parameters["num_features"].toInt();
		int clip_size=CLP.named_parameters["clip_size"].toInt();
		float ks_threshold=CLP.named_parameters["ks_threshold"].toFloat();
		int K_init=CLP.named_parameters["K_init"].toInt();

		if ((input_path.isEmpty())||(firings_path.isEmpty())||(output_path.isEmpty())) {cluster_usage(); return -1;}
		if (num_features==0) {cluster_usage(); return -1;}
		if (clip_size==0) {cluster_usage(); return -1;}
		if (ks_threshold==0) ks_threshold=1.4;
		if (K_init==0) K_init=25;

		if (!split_firings(input_path.toLatin1().data(),firings_path.toLatin1().data(),output_path.toLatin1().data(),num_features,clip_size,ks_threshold,K_init)) {
			printf("Error in cluster.\n");
			return -1;
		}
	}
    else if (command=="templates") {
        QString input_path=CLP.named_parameters["input"];
		QString firings_path=CLP.named_parameters["firings"];
        QString output_path=CLP.named_parameters["output"];
        int clip_size=CLP.named_parameters["clip_size"].toInt();

		if ((input_path.isEmpty())||(firings_path.isEmpty())||(output_path.isEmpty())) {templates_usage(); return -1;}
        if (clip_size==0) {templates_usage(); return -1;}

		if (!templates(input_path.toLatin1().data(),firings_path.toLatin1().data(),output_path.toLatin1().data(),clip_size)) {
            printf("Error in templates.\n");
            return -1;
        }
    }
    else if (command=="consolidate") {
		QString firings_path=CLP.named_parameters["firings"];
        QString templates_path=CLP.named_parameters["templates"];
        QString cluster_out_path=CLP.named_parameters["cluster_out"];
        QString templates_out_path=CLP.named_parameters["templates_out"];
        QString load_channels_out_path=CLP.named_parameters["load_channels_out"];
		float coincidence_threshold=CLP.named_parameters["coincidence_threshold"].toFloat();

		if ((firings_path.isEmpty())||(templates_path.isEmpty())) {consolidate_usage(); return -1;}
        if ((cluster_out_path.isEmpty())||(templates_out_path.isEmpty())) {consolidate_usage(); return -1;}
        if (load_channels_out_path.isEmpty()) {consolidate_usage(); return -1;}
		if (coincidence_threshold==0) {consolidate_usage(); return -1;}

		if (!consolidate(firings_path.toLatin1().data(),templates_path.toLatin1().data(),cluster_out_path.toLatin1().data(),templates_out_path.toLatin1().data(),load_channels_out_path.toLatin1().data(),coincidence_threshold)) {
            printf("Error in consolidate.\n");
            return -1;
        }
    }
    else if (command=="fit") {
        QString input_path=CLP.named_parameters["input"];
		QString firings_path=CLP.named_parameters["firings"];
        QString templates_path=CLP.named_parameters["templates"];
        QString cluster_out_path=CLP.named_parameters["cluster_out"];

		if ((input_path.isEmpty())||(firings_path.isEmpty())||(templates_path.isEmpty())) {fit_usage(); return -1;}
        if ((cluster_out_path.isEmpty())) {fit_usage(); return -1;}

		if (!fit(input_path.toLatin1().data(),templates_path.toLatin1().data(),firings_path.toLatin1().data(),cluster_out_path.toLatin1().data())) {
            printf("Error in fit.\n");
            return -1;
        }
    }
	else if (command=="extract_clips") {
		QString input_path=CLP.named_parameters["input"];
        QString detect_path=CLP.named_parameters["detect"];
        QString output_clips_path=CLP.named_parameters["output_clips"];
		int clip_size=CLP.named_parameters["clip_size"].toInt();

        if ((input_path.isEmpty())||(detect_path.isEmpty())) {extract_clips_usage(); return -1;}
        if (output_clips_path.isEmpty()) {extract_clips_usage(); return -1;}

        if (!extract_clips(input_path.toLatin1().data(),detect_path.toLatin1().data(),output_clips_path.toLatin1().data(),clip_size)) {
			printf("Error in extract_clips.\n");
			return -1;
		}
	}
    else if (command=="extract_channels") {
        QString input_path=CLP.named_parameters["input"];
        QString output_path=CLP.named_parameters["output"];
        QList<QString> channels_str=CLP.named_parameters["channels"].split(",");
        QList<int> channels;
        for (int i=0; i<channels_str.count(); i++) {
            channels << (channels_str[i].trimmed()).toInt();
        }

        if (input_path.isEmpty()) {extract_channels_usage(); return -1;}
        if (output_path.isEmpty()) {extract_channels_usage(); return -1;}
        if (channels.isEmpty()) {extract_channels_usage(); return -1;}

        if (!extract_channels(input_path.toLatin1().data(),output_path.toLatin1().data(),channels)) {
            printf("Error in extract_channels.\n");
            return -1;
        }
    }
    else if (command=="remove_artifacts") {
        QString input_path=CLP.named_parameters["input"];
        QString output_path=CLP.named_parameters["output"];
        float threshold=CLP.named_parameters["threshold"].toFloat();
        int normalize=CLP.named_parameters.value("normalize","1").toInt();
        int exclude_interval=CLP.named_parameters["exclude_interval"].toInt();

        if (input_path.isEmpty()) {printf("input path is empty.\n"); remove_artifacts_usage(); return -1;}
        if (output_path.isEmpty()) {printf("output path is empty.\n"); remove_artifacts_usage(); return -1;}
        if (threshold==0) {printf("threshold is zero.\n"); remove_artifacts_usage(); return -1;}
        if (exclude_interval==0) {printf("exclude_interval is zero.\n"); remove_artifacts_usage(); return -1;}

        if (!remove_artifacts(input_path.toLatin1().data(),output_path.toLatin1().data(),threshold,normalize,exclude_interval)) {
            printf("Error in remove_artifacts.\n");
            return -1;
        }
    }
    else if (command=="create_clips_file") {
        QString input_path=CLP.named_parameters["input"];
		QString firings_path=CLP.named_parameters["firings"];
        QString output_path=CLP.named_parameters["output"];
        QString index_out_path=CLP.named_parameters["index_out"];
        int clip_size=CLP.named_parameters["clip_size"].toInt();

		if ((input_path.isEmpty())||(firings_path.isEmpty())) {create_clips_file_usage(); return -1;}
        if ((output_path.isEmpty())||(index_out_path.isEmpty())) {create_clips_file_usage(); return -1;}

		if (!create_clips_file(input_path.toLatin1().data(),firings_path.toLatin1().data(),output_path.toLatin1().data(),index_out_path.toLatin1().data(),clip_size)) {
            printf("Error in create_clips_file.\n");
            return -1;
        }
    }
	else if (command=="assemble_firings_file") {
        QString input_detect_path=CLP.named_parameters["input_detect"];
        QString input_labels_path=CLP.named_parameters["input_labels"];
		QString output_firings_path=CLP.named_parameters["output_firings"];

		if ((input_detect_path.isEmpty())||(input_labels_path.isEmpty())) {assemble_firings_file_usage(); return -1;}
		if (output_firings_path.isEmpty()) {assemble_firings_file_usage(); return -1;}

		if (!assemble_firings_file(input_detect_path.toLatin1().data(),input_labels_path.toLatin1().data(),output_firings_path.toLatin1().data())) {
			printf("Error in assemble_firings_file.\n");
            return -1;
        }
    }
	else if (command=="cross_correlograms") {
		QString firings_path=CLP.named_parameters["firings"];
		QString output_path=CLP.named_parameters["output"];
		int max_dt=CLP.named_parameters["max_dt"].toInt();

		if ((firings_path.isEmpty())||(output_path.isEmpty())) {cross_correlograms_usage(); return -1;}
		if (max_dt==0) {cross_correlograms_usage(); return -1;}

		if (!cross_correlograms(firings_path.toLatin1(),output_path.toLatin1(),max_dt)) {
			printf("Error in cross_correlograms.\n");
			return -1;
		}
	}
	else if (command=="confusion_matrix") {
		QString firings1_path=CLP.named_parameters["firings1"];
		QString firings2_path=CLP.named_parameters["firings2"];
		QString output_path=CLP.named_parameters["output"];
		int max_matching_offset=CLP.named_parameters["max_matching_offset"].toInt();

		if ((firings1_path.isEmpty())||(firings2_path.isEmpty())||(output_path.isEmpty())) {confusion_matrix_usage(); return -1;}
		if (max_matching_offset==0) {confusion_matrix_usage(); return -1;}

		if (!confusion_matrix(firings1_path.toLatin1(),firings2_path.toLatin1(),output_path.toLatin1(),max_matching_offset)) {
			printf("Error in confusion_matrix.\n");
			return -1;
		}
	}
    else if (command=="branch_cluster_v1") {
        QString raw_path=CLP.named_parameters["raw"];
        QString detect_path=CLP.named_parameters["detect"];
        QString adjacency_matrix_path=CLP.named_parameters["adjacency_matrix"];
        QString firings_path=CLP.named_parameters["firings"];
        int clip_size=CLP.named_parameters.value("clip_size","0").toInt();
        int min_section_count=CLP.named_parameters.value("min_section_count","0").toInt();
        double section_increment=CLP.named_parameters.value("section_increment","0.5").toDouble();
        int num_features=CLP.named_parameters.value("num_features","0").toInt();

        if ((raw_path.isEmpty())||(detect_path.isEmpty())||(firings_path.isEmpty())) {branch_cluster_v1_usage(); return -1;}
        if (clip_size==0) {branch_cluster_v1_usage(); return -1;}
        if (min_section_count==0) {branch_cluster_v1_usage(); return -1;}
        if (section_increment==0) {branch_cluster_v1_usage(); return -1;}
        if (num_features==0) {branch_cluster_v1_usage(); return -1;}

        Branch_Cluster_Opts opts;
        opts.clip_size=clip_size;
        opts.min_section_count=min_section_count;
        opts.num_features=num_features;
        opts.section_increment=section_increment;
        if (!branch_cluster_v1(raw_path.toLatin1().data(),detect_path.toLatin1().data(),adjacency_matrix_path.toLatin1().data(),firings_path.toLatin1().data(),opts)) {
            printf("Error in branch_cluster_v1.\n");
            return -1;
        }
    }
	else if (command=="copy") {
		QString input_path=CLP.named_parameters["input"];
		QString output_path=CLP.named_parameters["output"];

		if ((input_path.isEmpty())||(output_path.isEmpty())) {copy_usage(); return -1;}
		if (!QFile::exists(input_path)) {
			printf("Error in copy... input file does not exist.\n");
			return -1;
		}
		if (QFile::exists(output_path)) QFile::remove(output_path);
		if (!QFile::copy(input_path,output_path)) {
			printf("Error in copy... unable to copy file.\n");
			return -1;
		}
	}
	else {
		printf("Unknown command: %s\n",command.toLatin1().data());
		return -1;
	}

	PT.reportProcessCompleted(CLP);

	printf("Elapsed time for %s: %.2f seconds\n",command.toLatin1().data(),timer.elapsed()*1.0/1000);

	return 0;
}
