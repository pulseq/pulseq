#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(pulseq, 1.1.0, kelvin.layton@uniklinik-freiburg.de)
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([foreign])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# C++ project
AC_LANG(C++)

# Checks for libraries.

# Checks for header files.

# Check for documentation option

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_INLINE
AC_ARG_ENABLE([documentation],
     [  --enable-documentation  turn on documentation generation (default: no)],
     [case "${enableval}" in
       yes) documentation=true ;;
       no)  documentation=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-documentation]) ;;
     esac],[documentation=false])
     AM_CONDITIONAL([DOCUMENTATION], [test x$documentation = xtrue])

if test x$documentation = xtrue; then

# Check for doxygen (required for documentation)
build_docs=yes
DOXYGEN=
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"; then
  AC_MSG_WARN([doxygen not found, not building code documentation])
  build_docs=no
fi

# Check for pdflatex
build_latex=yes
PDFLATEX=
AC_CHECK_PROGS([PDFLATEX],[pdflatex])
if test -z "$PDFLATEX"; then
  AC_MSG_WARN([pdflatex not found, not building specification document])
  build_latex=no
fi

fi
AM_CONDITIONAL(BUILD_DOCS, test x$build_docs = xyes)
AM_CONDITIONAL(BUILD_LATEX, test x$build_latex = xyes)

# Check for python (required for testing with "make check")
build_tests=yes
PYTHON=
AC_CHECK_PROGS([PYTHON], [python])
if test -z "$PYTHON"; then
  AC_MSG_WARN([python not found, cannot run test suite])
  build_tests=no
fi
AM_CONDITIONAL(BUILD_TESTS, test x$build_tests = xyes)

# Check for matlab
build_matlab_tests=yes
MATLAB=
AC_CHECK_PROGS([MATLAB], [matlab])
if test -z "$MATLAB"; then
  AC_MSG_WARN([matlab not found, cannot run matlab tests])
  build_matlab_tests=no
fi
AM_CONDITIONAL(BUILD_MATLAB_TESTS, test x$build_matlab_tests = xyes)

# Define output files
AC_CONFIG_FILES([Makefile src/Makefile matlab/Makefile doc/Makefile doc/doxy.conf])

AC_OUTPUT
