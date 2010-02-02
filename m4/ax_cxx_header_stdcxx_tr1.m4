# ===========================================================================
#    http://www.nongnu.org/autoconf-archive/ax_cxx_header_stdcxx_tr1.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_HEADER_STDCXX_TR1
#
# DESCRIPTION
#
#   Check for library coverage of the TR1 standard.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.
##
##  Modified by A. Gordon (gordon@cshl.edu), 1-Feb-2010
##  Removed unused header files (which can't be found on some TR1 gcc-4.2.4 on CentOS 5.4)
##

#serial 5

AU_ALIAS([AC_CXX_HEADER_STDCXX_TR1], [AX_CXX_HEADER_STDCXX_TR1])
AC_DEFUN([AX_CXX_HEADER_STDCXX_TR1], [
  AC_CACHE_CHECK(for ISO C++ TR1 include files,
  ax_cv_cxx_stdcxx_tr1,
  [AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  AC_TRY_COMPILE([
  #include <tr1/unordered_map>
  ],,
  ax_cv_cxx_stdcxx_tr1=yes, ax_cv_cxx_stdcxx_tr1=no)
  AC_LANG_RESTORE
  ])
  if test "$ax_cv_cxx_stdcxx_tr1" = yes; then
    AC_DEFINE(STDCXX_TR1_HEADERS,,[Define if ISO C++ TR1 header files are present. ])
  fi
])
