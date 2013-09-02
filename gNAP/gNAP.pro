#-------------------------------------------------
#
# Project created by QtCreator 2013-09-02T10:04:16
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = gNAP
TEMPLATE = app
RC_ICONS =icon.ico


SOURCES += main.cpp\
        widget.cpp \
    console.cpp \
    result.cpp \
    Code/Sequence.cpp \
    Code/SBOL.cpp \
    Code/RandSeq.cpp \
    Code/PSO.cpp \
    Code/ModleNetwork.cpp \
    Code/GRN.cpp \
    Code/GetReady.cpp \
    Code/GeneIM.cpp

HEADERS  += widget.h \
    console.h \
    result.h \
    Code/Sequence.h \
    Code/SBOL.h \
    Code/RandSeq.h \
    Code/PSO.h \
    Code/ModleNetwork.h \
    Code/GRN.h \
    Code/GetReady.h \
    Code/GeneIM.h \
    Code/define.h

FORMS    += widget.ui \
    console.ui \
    result.ui

RESOURCES += \
    gNAP.qrc
