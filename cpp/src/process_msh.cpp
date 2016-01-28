#include "process_msh.h"
#include "textfile.h"
#include <QCoreApplication>
#include <QFileInfo>
#include <QProcess>

int process_msh(const QString &path,int argc,char *argv[])
{
    QString txt=read_text_file(path);

    QString check_exit_status="rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi\n";

    QString txt2="";
    txt2+=QString("\n\nmountainview=%1/../mountainview/bin/mountainview\n\n").arg(qApp->applicationDirPath());
    txt2+=QString("\n\nscriptdir=%1\n\n").arg(QFileInfo(path).absolutePath());
    QList<QString> lines=txt.split("\n");
    QString buffer="";
    bool in_mscmd=false;
    for (int i=0; i<lines.count(); i++) {
        QString line=lines[i];
        if (line.indexOf("mscmd ")==0) {
            if (in_mscmd) {
                txt2+="\n";
                txt2+=check_exit_status;
                txt2+="echo\n";
            }
            txt2+=QString("%1/mountainsort ").arg(qApp->applicationDirPath())+line.mid(QString("mscmd ").count())+" ";
            in_mscmd=true;
        }
        else if (!line.trimmed().startsWith("--")) {
            if (in_mscmd) {
                txt2+="\n";
                txt2+=check_exit_status;
                txt2+="echo\n";
            }
            in_mscmd=false;
            txt2+=line+"\n";
        }
        else {
            txt2+=line;
            if (!in_mscmd) txt2+="\n";
        }
    }
    if (in_mscmd) {
        txt2+="\n";
        txt2+=check_exit_status;
        txt2+="echo\n";
    }

    write_text_file(path+".sh",txt2);

    QStringList args; args << path+".sh";
    for (int i=2; i<argc; i++) {
        args << argv[i];
    }
    QProcess::execute("/bin/bash",args);

    return 0;
}
