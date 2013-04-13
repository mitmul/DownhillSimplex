TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    downhill.cpp

# Eigen
INCLUDEPATH += /usr/local/include/eigen3

HEADERS += \
    downhill.h
