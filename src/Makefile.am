## Process this file with automake to produce Makefile.in

## Created by Anjuta

AM_CPPFLAGS = \
	-DPACKAGE_LOCALE_DIR=\""$(localedir)"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir)"\" \
	-DPACKAGE_DATA_DIR=\""$(pkgdatadir)"\" \
	$(GNU_R_SPHERE_CFLAGS)

AM_CFLAGS =\
	 -Wall\
	 -g

bin_PROGRAMS = gnu_r_sphere

gnu_r_sphere_SOURCES = \
	main.cc

gnu_r_sphere_LDFLAGS = 

gnu_r_sphere_LDADD = $(GNU_R_SPHERE_LIBS)


if NATIVE_WIN32
gnu_r_sphere_LDFLAGS += -mwindows
endif
