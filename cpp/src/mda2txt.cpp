#include "mda2txt.h"
#include "diskreadmda.h"
#include "textfile.h"
#include <stdio.h>

bool mda2txt(const char *mda_file, const char *txt_file, const mda2txt_opts &opts)
{
    DiskReadMda X; X.setPath(mda_file);

    QString txt;
    if (opts.transpose) {
        if (X.N1()>opts.max_cols) {
            printf("Too many columns!\n");
            return false;
        }
        if (X.N2()>opts.max_rows) {
            printf("Too many rows!\n");
            return false;
        }
        for (int i=0; i<X.N2(); i++) {
            QString line;
            for (int j=0; j<X.N1(); j++) {
                line+=QString("%1").arg(X.value(j,i));
                if (j+1<X.N1()) line+=QString("%1").arg(opts.delimiter);
            }
            txt+=line+"\n";
        }
    }
    else {
        if (X.N1()>opts.max_rows) {
            printf("Too many rows!\n");
            return false;
        }
        if (X.N2()>opts.max_cols) {
            printf("Too many columns!\n");
            return false;
        }
        for (int i=0; i<X.N1(); i++) {
            QString line;
            for (int j=0; j<X.N2(); j++) {
                line+=QString("%1").arg(X.value(i,j));
                if (j+1<X.N1()) line+=QString("%1").arg(opts.delimiter);
            }
            txt+=line+"\n";
        }
    }
    write_text_file(txt,txt_file);
    return true;
}
