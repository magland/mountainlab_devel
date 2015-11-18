QT += core
QT -= gui

CONFIG -= app_bundle #Please apple, don't make a bundle today

DESTDIR = ../../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainsort
TEMPLATE = app

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
    cluster.h
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
    cluster.cpp

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
QMAKE_CXXFLAGS += -fopenmp -std=c++11
LIBS += -fopenmp -lfftw3 -lfftw3_threads
