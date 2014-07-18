AC_INIT(parsnp,1.0)
AM_INIT_AUTOMAKE(parsnp,1.0)
AC_PROG_CC(gcc)

AC_ARG_WITH(libmuscle, [  --with-libmuscle=<path/to/libmuscle>     libMUSCLE install dir (default: `pwd`/muscle)])

if test "$with_libmuscle" == ""
then
	with_libmuscle=`pwd`/muscle
fi

AC_LANG(C++)

CPPFLAGS="-I$with_libmuscle"
AC_CHECK_HEADER(libMUSCLE/muscle.h, [result=1], [result=0])

if test $result == 0
then
	AC_MSG_ERROR([libMUSCLE  headers not found.])
fi

AC_SUBST(libmuscle, $with_libmuscle)
AC_OUTPUT(Makefile src/Makefile)