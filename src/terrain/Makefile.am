include $(top_srcdir)/Makefile.common

AM_CFLAGS += $(WARN_CFLAGS)
AM_CPPFLAGS += $(WARN_CPPFLAGS)
AM_CXXFLAGS += $(STD_CXXFLAGS) $(WARN_CXXFLAGS)

AM_CPPFLAGS += -I$(srcdir)/..

noinst_LIBRARIES = libterrain.a
noinst_HEADERS = \
	Terrain.h \
	TerrainNode.h \
	TerrainNoise.h

libterrain_a_SOURCES = \
	Terrain.cpp \
	TerrainNode.cpp \
	TerrainHeightJSON.cpp \
	TerrainColourJSON.cpp

AM_CPPFLAGS += -isystem @top_srcdir@/contrib
if !HAVE_LUA
AM_CPPFLAGS += -isystem @top_srcdir@/contrib/lua
endif
