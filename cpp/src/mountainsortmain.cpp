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
#include "split_clusters.h"
#include "templates.h"
#include "consolidate.h"
#include "extract_clips.h"
#include "get_principal_components.h"
#include "fit.h"
#include "cross_correlograms.h"
#include "confusion_matrix.h"
#include <QTime>

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
        P.version="0.13";
		PT.registerProcessor(P);
	}
    {
        PTProcessor P;
        P.command="features";
        P.input_file_pnames << "input";
        P.input_file_pnames << "detect";
        P.input_file_pnames << "adjacency";
        P.output_file_pnames << "output";
        P.version="0.17";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="cluster";
        P.input_file_pnames << "input";
        P.output_file_pnames << "output";
        P.version="0.16";
        PT.registerProcessor(P);
    }
	{
		PTProcessor P;
		P.command="split_clusters";
		P.input_file_pnames << "input";
		P.input_file_pnames << "cluster";
		P.output_file_pnames << "output";
        P.version="0.12";
		PT.registerProcessor(P);
	}
    {
        PTProcessor P;
        P.command="templates";
        P.input_file_pnames << "input";
        P.input_file_pnames << "cluster";
        P.output_file_pnames << "output";
		P.version="0.16";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="consolidate";
        P.input_file_pnames << "cluster";
        P.input_file_pnames << "templates";
        P.output_file_pnames << "cluster_out";
        P.output_file_pnames << "templates_out";
        P.output_file_pnames << "load_channels_out";
		P.version="0.39";
        PT.registerProcessor(P);
    }
    {
        PTProcessor P;
        P.command="fit";
        P.input_file_pnames << "input";
        P.input_file_pnames << "cluster";
        P.output_file_pnames << "templates";
        P.output_file_pnames << "cluster_out";
        P.version="0.13";
        PT.registerProcessor(P);
    }
	{
		PTProcessor P;
		P.command="extract_clips";
		P.input_file_pnames << "input";
		P.input_file_pnames << "cluster";
		P.output_file_pnames << "output";
		P.output_file_pnames << "index_out";
        P.version="0.11";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="cross_correlograms";
		P.input_file_pnames << "clusters";
		P.output_file_pnames << "output";
		P.version="0.10";
		PT.registerProcessor(P);
	}
	{
		PTProcessor P;
		P.command="confusion_matrix";
		P.input_file_pnames << "clusters1" << "clusters2";
		P.output_file_pnames << "output";
		P.version="0.11";
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
    printf("mountainsort features --input=in.mda --detect=detect.mda --adjacency=adjacency.mda --output=out.mda --num_features=6 --clip_size=100\n");
}

void cluster_usage() {
    printf("mountainsort cluster --input=in.mda --output=out.mda\n");
}

void split_clusters_usage() {
	printf("mountainsort split_clusters --input=in.mda --cluster=cluster.mda --output=out.mda --num_features=3 --clip_size=100\n");
}

void templates_usage() {
    printf("mountainsort templates --input=in.mda --cluster=cluster.mda --output=out.mda --clip_size=100\n");
}

void consolidate_usage() {
	printf("mountainsort consolidate --cluster=cluster.mda --templates=templates.mda --cluster_out=cluster_out.mda --templates_out=templates_out.mda --load_channels_out=load_channels.mda\n");
}

void fit_usage() {
    printf("mountainsort fit --input=input.mda --cluster=cluster.mda --templates=templates.mda --cluster_out=cluster_out.mda\n");
}

void extract_clips_usage() {
	printf("mountainsort extract_clips --input=raw.mda --cluster=cluster.mda --output=clips.mda --index_out=clips_index.mda --clip_size=100\n");
}

void cross_correlograms_usage() {
	printf("mountainsort cross_correlograms --clusters=clusters.mda --output=cross_correlograms.mda --max_dt=1500\n");
}

void confusion_matrix_usage() {
	printf("mountainsort confusion_matrix --clusters1=clusters1.mda --clusters2=clusters2.mda --output=confusion_matrix.mda --max_matching_offset=3\n");
}

/*

#include "armadillo"
void test_armadillo() {
    std::complex<double> x=0.1;
    std::complex<double> y=arma::arma_acos(x);
    printf("test armadillo: %g + %gi\n",y.real(),y.imag());
}

#include "pca.h"
void test_pca() {

    //float feps=std::numeric_limits<float>::epsilon();

    const int nvar = 4;
    stats::pca pca(nvar);

    const std::vector<double> record1 = {1, 2.5, 42, 7};
    const std::vector<double> record2 = {3, 4.2, 90, 7};
    const std::vector<double> record3 = {456, 444, 0, 7};
    pca.add_record(record1);
    pca.add_record(record2);
    pca.add_record(record3);

    pca.solve();

    //const float factor = 10;
    const std::vector<double> prinvec1 = pca.get_principal(0);
    const std::vector<double> exp_prinvec1 = {-2.10846198e+02, -2.13231575e+02, 4.24077773e+02};
    //assert_approx_equal_containers(exp_prinvec1, prinvec1, feps*factor, SPOT);

    const std::vector<double> prinvec2 = pca.get_principal(1);
    const std::vector<double> exp_prinvec2 = {-2.40512596e+01, 2.39612385e+01, 9.00211615e-02};
    //assert_approx_equal_containers(exp_prinvec2, prinvec2, feps*factor, SPOT);

    const std::vector<double> prinvec3 = pca.get_principal(2);
    const std::vector<double> exp_prinvec3 = {0, 0, 0};
    //assert_approx_equal_containers(exp_prinvec3, prinvec3, feps*factor, SPOT);

    const std::vector<double> prinvec4 = pca.get_principal(3);
    const std::vector<double> exp_prinvec4 = {0, 0, 0};
    //assert_approx_equal_containers(exp_prinvec4, prinvec4, feps*factor, SPOT);

    qDebug() << prinvec1[0] << prinvec1[1] << prinvec1[2];
    qDebug() << exp_prinvec1[0] << exp_prinvec1[1] << exp_prinvec1[2];
}
void test_pca_2() {

    //float feps=std::numeric_limits<float>::epsilon();

    const int nvar = 4;
    stats::pca pca(nvar);

    const std::vector<double> record1 = {1, 2,3,4.0};
    const std::vector<double> record2 = {2,4,6,8.8};
    const std::vector<double> record3 = {3,6,9,12.4};
    const std::vector<double> record4 = {1, 2,3,4.0};
    const std::vector<double> record5 = {2,4,6,8.8};
    const std::vector<double> record6 = {3,6,9,12.4};
    pca.add_record(record1);
    pca.add_record(record2);
    pca.add_record(record3);
    pca.add_record(record4);
    pca.add_record(record5);
    pca.add_record(record6);

    pca.solve();

    const std::vector<double> prinvec1 = pca.get_principal(0);
    const std::vector<double> prinvec2 = pca.get_principal(1);
    const std::vector<double> prinvec3 = pca.get_principal(2);
    const std::vector<double> prinvec4 = pca.get_principal(3);

    qDebug() << prinvec1[0] << prinvec1[1] << prinvec1[2] << prinvec1[3] << prinvec1[4] << prinvec1[5];
    qDebug() << prinvec2[0] << prinvec2[1] << prinvec2[2] << prinvec2[3] << prinvec2[4] << prinvec2[5];
    qDebug() << prinvec3[0] << prinvec3[1] << prinvec3[2] << prinvec3[3] << prinvec3[4] << prinvec3[5];
    qDebug() << prinvec4[0] << prinvec4[1] << prinvec4[2] << prinvec4[3] << prinvec4[4] << prinvec4[5];


	int ncomp=4;
	double *features=(double *)malloc(sizeof(double)*ncomp*6);
	double *data=(double *)malloc(sizeof(double)*nvar*6);
	int ii=0;
	for (int jj=0; jj<nvar; jj++) data[ii+jj]=record1[jj];
	ii+=nvar;
	for (int jj=0; jj<nvar; jj++) data[ii+jj]=record2[jj];
	ii+=nvar;
	for (int jj=0; jj<nvar; jj++) data[ii+jj]=record3[jj];
	ii+=nvar;
	for (int jj=0; jj<nvar; jj++) data[ii+jj]=record4[jj];
	ii+=nvar;
	for (int jj=0; jj<nvar; jj++) data[ii+jj]=record5[jj];
	ii+=nvar;
	for (int jj=0; jj<nvar; jj++) data[ii+jj]=record6[jj];
	ii+=nvar;
	get_pca_features(nvar,6,ncomp,features,data);
	printf("@@@@@@@@@@@@@@@@@@@@@@\n");
	for (int n=0; n<6; n++) {
		for (int m=0; m<nvar; m++) {
			printf("%g, ",data[m+nvar*n]);
		}
		printf("\n");
	}
	printf("@@@@@@@@@@@@@@@@@@@@@@\n");
	for (int k=0; k<ncomp; k++) {
		for (int n=0; n<6; n++) {
			printf("%g, ",features[k+ncomp*n]);
		}
		printf("\n");
	}
	free(features);
	free(data);
}
*/

int main(int argc,char *argv[]) {

	QCoreApplication app(argc,argv); //important for qApp->applicationDirPath() in processtracker

//	test_armadillo();
//	test_pca_2();
//	return 0;

	CLParams CLP;
	QStringList required;
	CLP=get_command_line_params(argc,argv,required);

	ProcessTracker PT;
	register_processors(PT);

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

        if (!cluster(input_path.toLatin1().data(),output_path.toLatin1().data())) {
            printf("Error in cluster.\n");
            return -1;
        }
    }
	else if (command=="split_clusters") {
		QString input_path=CLP.named_parameters["input"];
		QString cluster_path=CLP.named_parameters["cluster"];
		QString output_path=CLP.named_parameters["output"];
		int num_features=CLP.named_parameters["num_features"].toInt();
		int clip_size=CLP.named_parameters["clip_size"].toInt();

		if ((input_path.isEmpty())||(cluster_path.isEmpty())||(output_path.isEmpty())) {cluster_usage(); return -1;}
		if (num_features==0) {cluster_usage(); return -1;}
		if (clip_size==0) {cluster_usage(); return -1;}

		if (!split_clusters(input_path.toLatin1().data(),cluster_path.toLatin1().data(),output_path.toLatin1().data(),num_features,clip_size)) {
			printf("Error in cluster.\n");
			return -1;
		}
	}
    else if (command=="templates") {
        QString input_path=CLP.named_parameters["input"];
        QString cluster_path=CLP.named_parameters["cluster"];
        QString output_path=CLP.named_parameters["output"];
        int clip_size=CLP.named_parameters["clip_size"].toInt();

        if ((input_path.isEmpty())||(output_path.isEmpty())) {templates_usage(); return -1;}
        if (clip_size==0) {templates_usage(); return -1;}

        if (!templates(input_path.toLatin1().data(),cluster_path.toLatin1().data(),output_path.toLatin1().data(),clip_size)) {
            printf("Error in templates.\n");
            return -1;
        }
    }
    else if (command=="consolidate") {
        QString cluster_path=CLP.named_parameters["cluster"];
        QString templates_path=CLP.named_parameters["templates"];
        QString cluster_out_path=CLP.named_parameters["cluster_out"];
        QString templates_out_path=CLP.named_parameters["templates_out"];
        QString load_channels_out_path=CLP.named_parameters["load_channels_out"];

        if ((cluster_path.isEmpty())||(templates_path.isEmpty())) {consolidate_usage(); return -1;}
        if ((cluster_out_path.isEmpty())||(templates_out_path.isEmpty())) {consolidate_usage(); return -1;}
        if (load_channels_out_path.isEmpty()) {consolidate_usage(); return -1;}

		if (!consolidate(cluster_path.toLatin1().data(),templates_path.toLatin1().data(),cluster_out_path.toLatin1().data(),templates_out_path.toLatin1().data(),load_channels_out_path.toLatin1().data())) {
            printf("Error in consolidate.\n");
            return -1;
        }
    }
    else if (command=="fit") {
        QString input_path=CLP.named_parameters["input"];
        QString cluster_path=CLP.named_parameters["cluster"];
        QString templates_path=CLP.named_parameters["templates"];
        QString cluster_out_path=CLP.named_parameters["cluster_out"];

        if ((input_path.isEmpty())||(cluster_path.isEmpty())||(templates_path.isEmpty())) {fit_usage(); return -1;}
        if ((cluster_out_path.isEmpty())) {fit_usage(); return -1;}

        if (!fit(input_path.toLatin1().data(),templates_path.toLatin1().data(),cluster_path.toLatin1().data(),cluster_out_path.toLatin1().data())) {
            printf("Error in fit.\n");
            return -1;
        }
    }
	else if (command=="extract_clips") {
		QString input_path=CLP.named_parameters["input"];
		QString cluster_path=CLP.named_parameters["cluster"];
		QString output_path=CLP.named_parameters["output"];
		QString index_out_path=CLP.named_parameters["index_out"];
		int clip_size=CLP.named_parameters["clip_size"].toInt();

		if ((input_path.isEmpty())||(cluster_path.isEmpty())) {extract_usage(); return -1;}
		if ((output_path.isEmpty())||(index_out_path.isEmpty())) {extract_usage(); return -1;}

		if (!extract_clips(input_path.toLatin1().data(),cluster_path.toLatin1().data(),output_path.toLatin1().data(),index_out_path.toLatin1().data(),clip_size)) {
			printf("Error in extract_clips.\n");
			return -1;
		}
	}
	else if (command=="cross_correlograms") {
		QString clusters_path=CLP.named_parameters["clusters"];
		QString output_path=CLP.named_parameters["output"];
		int max_dt=CLP.named_parameters["max_dt"].toInt();

		if ((clusters_path.isEmpty())||(output_path.isEmpty())) {cross_correlograms_usage(); return -1;}
		if (max_dt==0) {cross_correlograms_usage(); return -1;}

		if (!cross_correlograms(clusters_path.toLatin1(),output_path.toLatin1(),max_dt)) {
			printf("Error in cross_correlograms.\n");
			return -1;
		}
	}
	else if (command=="confusion_matrix") {
		QString clusters1_path=CLP.named_parameters["clusters1"];
		QString clusters2_path=CLP.named_parameters["clusters2"];
		QString output_path=CLP.named_parameters["output"];
		int max_matching_offset=CLP.named_parameters["max_matching_offset"].toInt();

		if ((clusters1_path.isEmpty())||(clusters2_path.isEmpty())||(output_path.isEmpty())) {confusion_matrix_usage(); return -1;}
		if (max_matching_offset==0) {confusion_matrix_usage(); return -1;}

		if (!confusion_matrix(clusters1_path.toLatin1(),clusters2_path.toLatin1(),output_path.toLatin1(),max_matching_offset)) {
			printf("Error in confusion_matrix.\n");
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
