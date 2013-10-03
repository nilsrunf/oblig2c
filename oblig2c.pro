#-------------------------------------------------
#
# Project created by QtCreator 2013-10-03T12:40:52
#
#-------------------------------------------------

QT       = core

QT       -= gui

TARGET = oblig2c
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp
LIBS+= -larmadillo -llapack -lblas
