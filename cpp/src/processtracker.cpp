#include "processtracker.h"
#include <QCryptographicHash>
#include <QCoreApplication>
#include <QDir>
#include "qjson.h"
#include <sys/stat.h>

struct ProcessRecord {
	PTProcessor processor;
	QStringList input_file_codes;
	QMap<QString,QString> input_parameters;
	QStringList output_file_codes;
};

class ProcessTrackerPrivate {
public:
	ProcessTracker *q;
	QList<PTProcessor> m_registered_processors;

	PTProcessor find_processor_by_command(const QString &command,const QString &version);
	PTProcessor find_processor_by_command_only(const QString &command);
	QStringList compute_process_code_string_list(const PTProcessor &P,const CLParams &CLP);
	ProcessRecord find_process_record(const QString &code);
	QString compute_file_code(const QString &path);
	bool compare_string_lists(const QStringList &list1,const QStringList &list2);
	QString compute_hash_of_string_list(const QStringList &list);
	QString compute_hash_of_string(const QString &str);
	QString get_working_path();
	void store_process_record(const QString &code,const ProcessRecord &PR,const QStringList &input_file_paths,const QStringList &output_file_paths);
};

ProcessTracker::ProcessTracker()
{
	d=new ProcessTrackerPrivate;
	d->q=this;
}

ProcessTracker::~ProcessTracker()
{
	delete d;
}

QString read_text_file(const QString &path) {
	QFile f(path);
	if (!f.open(QFile::ReadOnly|QFile::Text)) return "";
	QTextStream in(&f);
	return in.readAll();
}
void write_text_file(const QString &path,const QString &txt) {
	QFile f(path);
	if (!f.open(QFile::WriteOnly|QFile::Text)) return;
	QTextStream out(&f);
	out << txt;
}

void ProcessTracker::registerProcessor(const PTProcessor &P)
{
	d->m_registered_processors << P;
}

bool ProcessTracker::processAlreadyCompleted(const CLParams &CLP)
{
	QString command=CLP.unnamed_parameters.value(0);
	PTProcessor P=d->find_processor_by_command_only(command);
	if (P.command.isEmpty()) {
		printf("Warning: processor not registered: %s\n",command.toLatin1().data());
		return false;
	}

	QStringList code_string_list=d->compute_process_code_string_list(P,CLP);
	if (code_string_list.isEmpty()) {
		return false;
	}
	QString code=d->compute_hash_of_string_list(code_string_list);

	ProcessRecord PR=d->find_process_record(code);
	if (PR.processor.command!=command) {
		return false;
	}

	if (PR.output_file_codes.count()!=P.output_file_pnames.count()) {
		return false;
	}
	QStringList codes0;
	for (int i=0; i<P.output_file_pnames.count(); i++) {
		QString pname=P.output_file_pnames[i];
		QString path=CLP.named_parameters.value(pname);
		codes0 << d->compute_file_code(path);
	}
	if (!d->compare_string_lists(PR.output_file_codes,codes0)) {
		return false;
	}
	return true;
}

void ProcessTracker::reportProcessCompleted(const CLParams &CLP)
{
	QString command=CLP.unnamed_parameters.value(0);
	PTProcessor P=d->find_processor_by_command_only(command);
	if (P.command.isEmpty()) {
		printf("Warning: unable to report process completed: %s\n",command.toLatin1().data());
		return;
	}
	ProcessRecord PR;
	PR.processor=P;
	QStringList keys=CLP.named_parameters.keys();
	foreach (QString key,keys) {
		if ((!P.input_file_pnames.contains(key))&&(!P.output_file_pnames.contains(key))) {
			PR.input_parameters[key]=CLP.named_parameters[key];
		}
	}
	QStringList input_file_paths;
	for (int i=0; i<P.input_file_pnames.count(); i++) {
		QString path0=CLP.named_parameters[P.input_file_pnames[i]];
		if (path0.isEmpty()) {
			printf("Warning: unable to report process completed: input file name is empty: %s.\n",P.input_file_pnames[i].toLatin1().data());
			return;
		}
		QString code0=d->compute_file_code(path0);
		if (code0.isEmpty()) {
			printf("Warning: unable to report process completed: input file code is empty: %s.\n",P.input_file_pnames[i].toLatin1().data());
			return;
		}
		input_file_paths << path0;
		PR.input_file_codes << code0;
	}
	QStringList output_file_paths;
	for (int i=0; i<P.output_file_pnames.count(); i++) {
		QString path0=CLP.named_parameters[P.output_file_pnames[i]];
		if (path0.isEmpty()) {
			printf("Warning: unable to report process completed: output file name is empty: %s.\n",P.output_file_pnames[i].toLatin1().data());
			return;
		}
		QString code0=d->compute_file_code(path0);
		if (code0.isEmpty()) {
			printf("Warning: unable to report process completed: output file code is empty: %s.\n",P.output_file_pnames[i].toLatin1().data());
			return;
		}
		output_file_paths << path0;
		PR.output_file_codes << code0;
	}

	QStringList code_string_list=d->compute_process_code_string_list(P,CLP);
	if (code_string_list.isEmpty()) {
		printf("Warning: unable to report process completed: code_string_list is empty.\n");
		return;
	}
	QString code=d->compute_hash_of_string_list(code_string_list);

	d->store_process_record(code,PR,input_file_paths,output_file_paths);
}

int ProcessTracker::processorCount()
{
	return d->m_registered_processors.count();
}

PTProcessor ProcessTracker::processor(int i)
{
	return d->m_registered_processors[i];
}

PTProcessor ProcessTracker::findProcessor(const QString &command)
{
	return d->find_processor_by_command_only(command);
}


PTProcessor ProcessTrackerPrivate::find_processor_by_command(const QString &command,const QString &version)
{
	for (int i=0; i<m_registered_processors.count(); i++) {
		if ((m_registered_processors[i].command==command)&&(m_registered_processors[i].version==version)) {
			return m_registered_processors[i];
		}
	}
	PTProcessor empty_processor;
	return empty_processor;
}

PTProcessor ProcessTrackerPrivate::find_processor_by_command_only(const QString &command)
{
	for (int i=0; i<m_registered_processors.count(); i++) {
		if (m_registered_processors[i].command==command) {
			return m_registered_processors[i];
		}
	}
	PTProcessor empty_processor;
	return empty_processor;
}

QStringList ProcessTrackerPrivate::compute_process_code_string_list(const PTProcessor &P, const CLParams &CLP)
{
	QStringList code_strings;
	code_strings << "--command--";
	code_strings << P.command;
	code_strings << "--version--";
	code_strings << P.version;
	code_strings << "--input files--";
	for (int i=0; i<P.input_file_pnames.count(); i++) {
		QString pname=P.input_file_pnames[i];
		QString code=compute_file_code(CLP.named_parameters.value(pname));
		if (code.isEmpty()) return QStringList();
		QString str=QString("%1:%2").arg(pname).arg(code);
		code_strings << str;
	}
	code_strings << "--parameters--";
	QStringList keys=CLP.named_parameters.keys();
	qSort(keys);
	foreach (QString key,keys) {
		if ((!P.input_file_pnames.contains(key))&&(!P.output_file_pnames.contains(key))) {
			QString str=QString("%1=%2").arg(key).arg(CLP.named_parameters.value(key));
			code_strings << str;
		}
	}
	return code_strings;
}

QStringList to_string_list(const QList<QVariant> &list) {
	QStringList ret;
	for (int i=0; i<list.count(); i++) {
		ret << list[i].toString();
	}
	return ret;
}

QList<QVariant> to_variant_list(const QStringList &list) {
	QList<QVariant> ret;
	for (int i=0; i<list.count(); i++) {
		ret << list[i];
	}
	return ret;
}

QMap<QString,QVariant> to_variant_map(const QMap<QString,QString> &map) {
	QMap<QString,QVariant> ret;
	QStringList keys=map.keys();
	foreach (QString key,keys) {
		ret[key]=map[key];
	}
	return ret;
}

QMap<QString,QString> to_string_map(const QMap<QString,QVariant> &map) {
	QMap<QString,QString> ret;
	QStringList keys=map.keys();
	foreach (QString key,keys) {
		ret[key]=map[key].toString();
	}
	return ret;
}

ProcessRecord ProcessTrackerPrivate::find_process_record(const QString &code)
{
	QString path0=get_working_path();
	ProcessRecord PR;
	QString path1=QString("%1/%2.process").arg(path0).arg(code);
	if (QFile::exists(path1)) {
		QString txt=read_text_file(path1);
		QVariant json=parseJSON(txt);
		QMap<QString,QVariant> X=json.toMap();
		QString command=X.value("command").toString();
		QString version=X.value("version").toString();
		PR.processor=find_processor_by_command(command,version);
		if (PR.processor.command.isEmpty()) return PR;
		if (PR.processor.command!=command) return PR;
		PR.input_file_codes=to_string_list(X.value("input_file_codes").toList());
		PR.input_parameters=to_string_map(X.value("input_parameters").toMap());
		PR.output_file_codes=to_string_list(X.value("output_file_codes").toList());
	}
	return PR;
}

void ProcessTrackerPrivate::store_process_record(const QString &code,const ProcessRecord &PR,const QStringList &input_file_paths,const QStringList &output_file_paths)
{
	QString path0=get_working_path();
	QString path1=QString("%1/%2.process").arg(path0).arg(code);
	QMap<QString,QVariant> X;
	X["command"]=PR.processor.command;
	X["version"]=PR.processor.version;
	X["input_file_codes"]=to_variant_list(PR.input_file_codes);
	X["input_parameters"]=to_variant_map(PR.input_parameters);
	X["output_file_codes"]=to_variant_list(PR.output_file_codes);
	X["input_file_paths"]=to_variant_list(input_file_paths);
	X["output_file_paths"]=to_variant_list(output_file_paths);
	QString txt=toJSON(X);
	write_text_file(path1,txt);
}


QString ProcessTrackerPrivate::compute_file_code(const QString &path)
{
	//the code comprises the device,inode,size, and modification time (in seconds)
	//note that it is not dependent on the file name
	struct stat SS;
	stat(path.toLatin1().data(),&SS);
	QString id_string=QString("%1:%2:%3:%4").arg(SS.st_dev).arg(SS.st_ino).arg(SS.st_size).arg(SS.st_mtim.tv_sec);
	return id_string;
}

bool ProcessTrackerPrivate::compare_string_lists(const QStringList &list1, const QStringList &list2)
{
	if (list1.count()!=list2.count()) return false;
	for (int i=0; i<list1.count(); i++) {
		if (list1[i]!=list2[i]) return false;
	}
	return true;
}

QString ProcessTrackerPrivate::compute_hash_of_string_list(const QStringList &list)
{
	QString str;
	str+=QString("&&&%1&&& ").arg(list.count());
	for (int i=0; i<list.count(); i++) {
		str+=QString("(%1)%2").arg(list[i].count()).arg(list[i]);
	}
	return compute_hash_of_string(str);
}

QString ProcessTrackerPrivate::compute_hash_of_string(const QString &str)
{
	QCryptographicHash hash(QCryptographicHash::Sha1);
	hash.addData(str.toLatin1());
	return QString(hash.result().toHex());
}

QString ProcessTrackerPrivate::get_working_path()
{
	QString path0=qApp->applicationDirPath();
	if (!QDir(path0).exists(".process_tracker")) {
		QDir(path0).mkdir(".process_tracker");
	}
	return path0+"/.process_tracker";
}

