QT += core gui
QT += script

CONFIG -= app_bundle #Please apple, don't make a bundle today

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets #We do want to support Qt5, but there is no reason not to use Qt4

DESTDIR = ../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainview
TEMPLATE = app

HEADERS += mountainviewwidget.h \
    histogramview.h \
    mvoverviewwidget.h \
    mvstatisticswidget.h \
    mvcrosscorrelogramswidget.h \
    mvunitwidget.h \
    diskarraymodelclipssubset.h
SOURCES += mountainviewmain.cpp mountainviewwidget.cpp \
    histogramview.cpp \
    mvoverviewwidget.cpp \
    mvstatisticswidget.cpp \
    mvcrosscorrelogramswidget.cpp \
    mvunitwidget.cpp \
    diskarraymodelclipssubset.cpp

HEADERS += get_command_line_params.h
SOURCES += get_command_line_params.cpp

INCLUDEPATH += spikespy/src
DEPENDPATH += spikespy/src #This DEPENDPATH is for Qt4
VPATH += spikespy/src #This VPATH is for Qt5

SOURCES += \
    sstimeseriesplot.cpp \
    plotarea.cpp \
    sstimeseriesview.cpp \
    sscontroller.cpp \
    sstimeserieswidget.cpp \
    mdaobject.cpp \
    diskarraymodel.cpp \
    sslabelsmodel1.cpp \
    mdaio.cpp \
    ssabstractview.cpp \
    sslabelview.cpp \
    ssabstractplot.cpp \
    sslabelplot.cpp \
    extractclipsdialog.cpp \
    cvcombowidget.cpp \
    diskreadmda.cpp \
    diskwritemda.cpp \
    memorymda.cpp \
    usagetracking.cpp \
    sscommon.cpp

HEADERS  += sstimeseriesplot.h \
    plotarea.h \
    sstimeseriesview.h \
    sstimeserieswidget.h \
    sscontroller.h \
    mdaobject.h \
    diskarraymodel.h \
    sscommon.h \
    mdaio.h \
    sslabelsmodel.h \
    sslabelsmodel1.h \
    ssabstractview.h \
    sslabelview.h \
    ssabstractplot.h \
    sslabelplot.h \
    extractclipsdialog.h \
    cvcombowidget.h \
    diskreadmda.h \
    diskwritemda.h \
    memorymda.h \
    usagetracking.h

SOURCES += cvwidget.cpp cvview.cpp affinetransformation.cpp cvcommon.cpp
HEADERS += cvwidget.h cvview.h affinetransformation.h cvcommon.h

HEADERS += mda.h textfile.h
SOURCES += mda.cpp textfile.cpp

FORMS += \
    extractclipsdialog.ui

INCLUDEPATH += firetrack/src
DEPENDPATH += firetrack/src
VPATH += firetrack/src
HEADERS += \
    ftelectrodearrayview.h \
    ftelectrodearraywidget.h \
    firetrackwidget.h \
    ftoptionswidget.h \
    fthelpwidget.h \
    ftplotoptions.h
SOURCES += \
    ftelectrodearrayview.cpp \
    ftelectrodearraywidget.cpp \
    firetrackwidget.cpp \
    ftoptionswidget.cpp \
    fthelpwidget.cpp \
    ftplotoptions.cpp

RESOURCES += mountainview.qrc
