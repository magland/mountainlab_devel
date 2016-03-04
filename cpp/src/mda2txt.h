#ifndef MDA2TXT
#define MDA2TXT

struct mda2txt_opts {
    char delimiter;
    bool transpose;
    long max_rows,max_cols;
};

bool mda2txt(const char *mda_file,const char *txt_file,const mda2txt_opts &opts);

#endif // MDA2TXT

