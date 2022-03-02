#
#    ESDM-PAV Analytical Kernels
#    Copyright (C) 2022 CMCC Foundation
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

AC_PREREQ([2.60])
AC_INIT([esdm-pav-analytics-kernels], [1.0.0], [ophidia-info@cmcc.it])

AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AM_PROG_CC_STDC
AM_PROG_CC_C_O
AC_PROG_CXX

LT_INIT([disable-static])

have_esdm=no
ESDM_CFLAGS=
ESDM_LIBS=
AC_ARG_WITH(esdm_path,
           [  --with-esdm-path=PATH  Set the esdm root dir],
           [
				have_esdm=yes
				ESDM_CFLAGS="-I$withval/include"
				ESDM_LIBS="-L$withval/lib -lesdm -Wl,-rpath -Wl,$withval/lib"
				AC_SUBST(ESDM_CFLAGS)
				AC_SUBST(ESDM_LIBS)
				AC_MSG_NOTICE([ESDM support enabled])
           ],
           [
				AC_MSG_ERROR([ESDM library not found! Aborting configure!])
           ]
           )
AM_CONDITIONAL([HAVE_ESDM], [test "x$have_esdm" = "xyes"])

OPT="-Wno-error -Wno-format-security"
case "${host}" in
        *-*-solaris*)   PLATFORM=SUN_OS
                        SAMPLE_EXTRA_LIBS="-lxnet -lsocket -lnsl"
                        ;;
        *-*-sysv5*)     PLATFORM=UNIXWARE ;;
        *-*-UnixWare*)  PLATFORM=UNIXWARE ;;
        *-*-unixware*)  PLATFORM=UNIXWARE ;;
        *-*-sco3.2v5*)  PLATFORM=OPENSERVER
                        SAMPLE_EXTRA_LIBS="-lsocket"
                        ;;
        *-*-linux*)     PLATFORM=LINUX
                                    if (test "x$debug" = "x1") || (test "x$debug" = "x2"); then
										OPT+=" -g3 -Wall -Wextra -O0"
									else
										OPT+=" -O3"
									fi
			            if test "x$optimization" = "xyes" ; then
			                 CFLAGS="-O3 -finline-functions -funswitch-loops -fgcse-after-reload -ftree-loop-linear -ftree-loop-im -ftree-loop-ivcanon -fivopts -ftree-vectorize -ftracer -fgcse-sm -fgcse-las -fgcse-after-reload -funroll-loops -freorder-blocks-and-partition -fweb"
			            fi
			            ;;
                *-*-cygwin*)    PLATFORM=CYGWIN
                        ;;
        *-*-freebsd*)   PLATFORM=FREEBSD ;;
        *-*-irix*)      PLATFORM=IRIX
                        ;;
        *-*-aix*)       PLATFORM=AIX
                      case "${host}" in
                         *-*-aix4.1*)   osver=AIX41;;
                         *-*-aix4.2*)   osver=AIX42;;
                         *-*-aix4.3*)   osver=AIX43;;
                         *-*-aix4.*)    osver=AIX43;;
                      esac;;
        *-sequent-*)    PLATFORM=PTX ;;
        *-*-hp*)        PLATFORM=HP_UX ;
                        case "${host}" in
                            *-*-hpux11*)  osver=HPUX11;;
                            *-*-hpux10*)  osver=HPUX10;;
                            *)            osver=HPUX11;;
                        esac;;
        *-*-mvs*)       PLATFORM=OS390 ;;
        *-*-os400*)     PLATFORM=OS400 ;;
        *-*-OS400*)     PLATFORM=OS400 ;;
        *-*-osf*)       PLATFORM=TRU64
				        OPT+=" -g3"
				        REENTRANT="1"
                               ;;
        *-apple-*)      PLATFORM=__APPLE__
				        OPT+=" -g3"
				        ;;
        *)              PLATFORM=UNKNOWN ;;
esac

AC_SUBST(SAMPLE_EXTRA_LIBS)
AC_SUBST(PLATFORM)
AC_SUBST(OPT)
AC_SUBST(REENTRANT)

AC_OUTPUT([Makefile src/Makefile])
