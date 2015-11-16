#include <stdlib.h>
#include <stdio.h>

#include "get_command_line_params.h"
#include "bandpass_filter.h"

int main(int argc,char *argv[]) {
	CLParams CLP;
	QStringList required;
	QStringList optional;
	optional << "input" << "output";
	optional << "samplefreq" << "freq_min" << "freq_max";
	CLP=get_command_line_params(argc,argv,required,optional);

	if (CLP.unnamed_parameters.count()!=1) {
		printf("No command specified.\n");
		return -1;
	}
	QString command=CLP.unnamed_parameters.value(0);
	QString input_path=CLP.named_parameters["input"];
	QString output_path=CLP.named_parameters["output"];

	if (command=="filter") {
		double samplefreq=CLP.named_parameters["samplefreq"].toDouble();
		double freq_min=CLP.named_parameters["freq_min"].toDouble();
		double freq_max=CLP.named_parameters["freq_max"].toDouble();

		bandpass_filter(input_path.toLatin1().data(),output_path.toLatin1().data(),samplefreq,freq_min,freq_max);
	}

	return 0;
}
