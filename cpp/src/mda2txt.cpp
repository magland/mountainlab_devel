#include "mda2txt.h"
#include "diskreadmda.h"
#include "textfile.h"
#include <stdio.h>

QString format_number(double num) {
    if (num==(long)num) {
        return QString::number(num,'f',0);
    }
    else {
        return QString::number(num,'f');
    }
}

bool mda2txt(const char *mda_file, const char *txt_file, const mda2txt_opts &opts)
{
    DiskReadMda X; X.setPath(mda_file);

    QString txt;
    if (opts.transpose) {
        if (X.N1()>opts.max_cols) {
            printf("Too many columns::: %ld>%ld\n",(long)X.N1(),opts.max_cols);
            return false;
        }
        if (X.N2()>opts.max_rows) {
            printf("Too many rows::: %ld>%ld\n",(long)X.N2(),opts.max_rows);
            return false;
        }
        for (int i=0; i<X.N2(); i++) {
            QString line;
            for (int j=0; j<X.N1(); j++) {
                line+=QString("%1").arg(format_number(X.value(j,i)));
                if (j+1<X.N1()) line+=QString("%1").arg(opts.delimiter);
            }
            txt+=line+"\n";
        }
    }
    else {
        if (X.N1()>opts.max_rows) {
            printf("Too many rows: %ld\n",(long)X.N1());
            return false;
        }
        if (X.N2()>opts.max_cols) {
            printf("Too many columns: %ld\n",(long)X.N2());
            return false;
        }
        for (int i=0; i<X.N1(); i++) {
            QString line;
            for (int j=0; j<X.N2(); j++) {
                line+=QString("%1").arg(format_number(X.value(i,j)));
                if (j+1<X.N1()) line+=QString("%1").arg(opts.delimiter);
            }
            txt+=line+"\n";
        }
    }
    printf("Writing %d bytes to %s...\n",txt.count(),txt_file);
    if (!write_text_file(txt_file,txt)) {
        printf("Unable to write file %s\n",txt_file);
        return false;
    }
    return true;
}
