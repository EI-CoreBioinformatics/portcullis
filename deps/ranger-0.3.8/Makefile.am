
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(top_srcdir)/ranger.pc

# Create portcullis library
lib_LTLIBRARIES = libranger.la

# Library inputs
libranger_la_LDFLAGS = -version-info 3:0:3
libranger_la_SOURCES =	\
	src/DataChar.cpp \
	src/Data.cpp \
	src/DataDouble.cpp \
	src/DataFloat.cpp \
	src/ForestClassification.cpp \
	src/Forest.cpp \
	src/ForestProbability.cpp \
	src/ForestRegression.cpp \
	src/ForestSurvival.cpp \
	src/TreeClassification.cpp \
	src/Tree.cpp \
	src/TreeProbability.cpp \
	src/TreeRegression.cpp \
	src/TreeSurvival.cpp \
	src/utility.cpp

library_includedir=$(includedir)/ranger-@PACKAGE_VERSION@/ranger
RI = include/ranger
library_include_HEADERS = \
	$(RI)/DataChar.h \
	$(RI)/DataDouble.h \
	$(RI)/DataFloat.h \
	$(RI)/Data.h \
	$(RI)/ForestClassification.h \
	$(RI)/Forest.h \
	$(RI)/ForestProbability.h \
	$(RI)/ForestRegression.h \
	$(RI)/ForestSurvival.h \
	$(RI)/globals.h \
	$(RI)/TreeClassification.h \
	$(RI)/Tree.h \
	$(RI)/TreeProbability.h \
	$(RI)/TreeRegression.h \
	$(RI)/TreeSurvival.h \
	$(RI)/utility.h \
	$(RI)/version.h

libranger_la_CPPFLAGS = \
	-isystem $(top_srcdir)/include \
	@AM_CPPFLAGS@

libranger_la_CXXFLAGS = -Wall -Wextra -Wno-unused-parameter @AM_CXXFLAGS@
