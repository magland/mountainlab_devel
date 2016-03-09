QT += core
QT -= gui

CONFIG -= app_bundle #Please apple, don't make a bundle today

DESTDIR = ../../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainsort
TEMPLATE = app

#DEFINES += USE_LAPACK
#LIBS += -llapack -llapacke

# we want to avoid using any 3rdparty libraries -- this is important! Because ease of installation is critical!
#LIBS += -larmadillo -lpca -L/usr/local/lib -Wl,-rpath,$(DEFAULT_LIB_INSTALL_PATH)

INCLUDEPATH += ../../mountainview/src/spikespy/src
DEPENDPATH += ../../mountainview/src/spikespy/src
VPATH += ../../mountainview/src/spikespy/src
HEADERS += diskreadmda.h \
    mda2txt.h \
    outlier_scores_v1.h
SOURCES += diskreadmda.cpp \
    mda2txt.cpp \
    outlier_scores_v1.cpp

HEADERS += \
    bandpass_filter.h \
    usagetracking.h \
    mdaio.h \
    processtracker.h \
    normalize_channels.h \
    whiten.h \
    pcasolver.h \
    extract.h \
    detect.h \
    features0.h \
    get_principal_components.h \
    cluster.h \
    templates.h \
    consolidate.h \
    split_firings.h \
    extract_clips.h \
    fit.h \
    cross_correlograms.h \
    confusion_matrix.h \
    isobranch.h \
    create_clips_file.h \
    process_msh.h \
    textfile.h \
    assemble_firings_file.h \
    extract_channels.h \
    remove_artifacts.h \
    branch_cluster_v1.h \
    remove_duplicates.h \
    get_sort_indices.h \
    remove_noise_subclusters.h \
    msutils.h

SOURCES += \
mountainsortmain.cpp \
    bandpass_filter.cpp \
    usagetracking.cpp \
    mdaio.cpp \
    processtracker.cpp \
    normalize_channels.cpp \
    whiten.cpp \
    pcasolver.cpp \
    extract.cpp \
    detect.cpp \
    features0.cpp \
    get_principal_components.cpp \
    cluster.cpp \
    templates.cpp \
    consolidate.cpp \
    split_firings.cpp \
    extract_clips.cpp \
    fit.cpp \
    cross_correlograms.cpp \
    confusion_matrix.cpp \
    isobranch.cpp \
    create_clips_file.cpp \
    process_msh.cpp \
    textfile.cpp \
    assemble_firings_file.cpp \
    extract_channels.cpp \
    remove_artifacts.cpp \
    branch_cluster_v1.cpp \
    remove_duplicates.cpp \
    get_sort_indices.cpp \
    remove_noise_subclusters.cpp \
    msutils.cpp

INCLUDEPATH += ../../isosplit
DEPENDPATH += ../../isosplit
VPATH += ../../isosplit
HEADERS += isocut.h mda.h jisotonic.h isosplit.h isosplit2.h
SOURCES += isocut.cpp mda.cpp jisotonic.cpp isosplit.cpp isosplit2.cpp

HEADERS += qjson.h
SOURCES += qjson.cpp
INCLUDEPATH += qjson
DEPENDPATH += qjson
VPATH += qjson
HEADERS += serializer.h serializerrunnable.h parser.h parserrunnable.h json_scanner.h json_parser.hh
SOURCES += serializer.cpp serializerrunnable.cpp parser.cpp parserrunnable.cpp json_scanner.cpp json_parser.cc

HEADERS += get_command_line_params.h
SOURCES += get_command_line_params.cpp

QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp
#-std=c++11   # AHB removed since not in GNU gcc 4.6.3

LIBS += -fopenmp -lfftw3 -lfftw3_threads
