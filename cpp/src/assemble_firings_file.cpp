#include "assemble_firings_file.h"

#include "mda.h"
#include <QList>

bool assemble_firings_file(const char *input_detect_path,const char *input_labels_path,const char *output_firings_path) {
    Mda D; D.read(input_detect_path);
    Mda L; L.read(input_labels_path);

    QList<int> inds;
    for (int i=0; i<D.N2(); i++) {
        if (L.value(0,i)>0) inds << i;
    }

    Mda C;
    C.allocate(3,inds.count());
    for (int i=0; i<inds.count(); i++) {
        C.setValue(D.value(0,inds[i]),0,i);
        C.setValue(D.value(1,inds[i]),1,i);
        C.setValue(L.value(0,inds[i]),2,i);
    }
	C.write(output_firings_path);
    return true;
}
