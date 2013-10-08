TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lconfig++

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
