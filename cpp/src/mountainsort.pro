QT += core
QT -= gui

CONFIG -= app_bundle #Please apple, don't make a bundle today

DESTDIR = ../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainsort
TEMPLATE = app

HEADERS += \ 
    bandpass_filter.h \
    usagetracking.h \
    mdaio.h
SOURCES += mountainsortmain.cpp \
    bandpass_filter.cpp \
    usagetracking.cpp \
    mdaio.cpp

HEADERS += get_command_line_params.h
SOURCES += get_command_line_params.cpp

QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -std=c++11
LIBS += -fopenmp -lfftw3 -lfftw3_threads
