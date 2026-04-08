lib_LTLIBRARIES = libbraiding.la
libbraiding_la_SOURCES = cbraid.cpp braiding.cpp braiding.h cbraid.h cbraid_implementation.h  cbraid_interface.h
libbraiding_la_LDFLAGS = -lm
include_HEADERS = braiding.h cbraid.h cbraid_implementation.h  cbraid_interface.h
AM_CPPFLAGS = -Wall -fPIC -I$(LIBCBRAID_INCLUDEDIR) $(CPPFLAGS_CLN)
# Use make DEBUG=1 to compile with debug info.
if CPPDEBUG
AM_CPPFLAGS += -g
else
# Only optimize when the DEBUG flag is unset.
AM_CPPFLAGS += -O
endif
