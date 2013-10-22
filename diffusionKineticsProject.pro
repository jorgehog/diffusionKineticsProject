TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lconfig++ -lpython2.7
INCLUDEPATH += /usr/include/python2.7
QMAKE_CXXFLAGS += -std=c++0x

SOURCES += main.cpp \
    diffusionscheme.cpp \
    expliciteuler.cpp \
    solver.cpp

HEADERS += \
    diffusionscheme.h \
    expliciteuler.h \
    Schemes.h \
    diffusionKinetics.h \
    solver.h \
    constants.h

OTHER_FILES += config.cfg

release {

    QMAKE_CXXFLAGS -= O2
    QMAKE_CXXFLAGS += -O3 -DARMA_NO_DEBUG

}
