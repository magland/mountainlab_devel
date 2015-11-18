#ifndef PROCESSTRACKER_H
#define PROCESSTRACKER_H

#include "get_command_line_params.h"

struct PTProcessor {
	QString command;
	QString version;
	QStringList input_file_pnames;
	QStringList output_file_pnames;
};

class ProcessTrackerPrivate;
class ProcessTracker
{
public:
	friend class ProcessTrackerPrivate;
	ProcessTracker();
	virtual ~ProcessTracker();
	void registerProcessor(const PTProcessor &P);
	bool processAlreadyCompleted(const CLParams &CLP);
	void reportProcessCompleted(const CLParams &CLP);
	int processorCount();
	PTProcessor processor(int i);
	PTProcessor findProcessor(const QString &command);

private:
	ProcessTrackerPrivate *d;
};

#endif // PROCESSTRACKER_H
