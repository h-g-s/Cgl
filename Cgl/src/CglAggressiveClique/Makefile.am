AUTOMAKE_OPTIONS = foreign

########################################################################
#                           libCglAggressiveClique                     #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglAggressiveClique.la

# List all source files for this library, including headers
libCglAggressiveClique_la_SOURCES = \
	CglAggressiveClique.cpp CglAggressiveClique.hpp

# This is for libtool (on Windows)
libCglAggressiveClique_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/..` $(COINUTILS_CFLAGS) $(OSI_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = CglAggressiveClique.hpp
