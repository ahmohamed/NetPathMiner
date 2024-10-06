AC_INIT([NetPathMiner]) dnl package name, version
AC_CONFIG_SRCDIR(DESCRIPTION) dnl make sure the dir is an R package

dnl A user-specifiable option

dnl Select an optional include path, from a configure option
dnl or from an environment variable.



: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP

AC_CHECK_LIB(xml2, xmlParseFile,, [XML=false])
AC_CHECK_LIB(sbml, readSBML, , [SBML=false])


AC_DEFUN([CHECK_XML],[
  AC_PATH_PROGS(XML_CONFIG, xml2-config)
  LIBXML_INCDIR="-DHAVE_XML `${XML_CONFIG} --cflags`"
  XML_LIBS=`${XML_CONFIG} --libs`
])

PKG_LIBS="${LIBS}"


if ${SBML}
  then
    echo 'libSBML found, finding configuration...';
    if type "pkg-config" > /dev/null
      then
          PKG_CONFIG="pkg-config"
          echo "pkg-config tool found"
      else
        AC_PATH_PROG(PKG_CONFIG, pkg-config)
    fi;

    SBMLPC="libsbml"
    if test -z `${PKG_CONFIG} --exists ${SBMLPC}`
      then
        SBMLPC=`locate libsbml.pc 2>/dev/null`
    fi;

    if test -n "${SBMLPC}"
      then
        echo "libsbml.pc found. Using pkg-config to extract it."
        SBML_INCDIR="-DHAVE_SBML `${PKG_CONFIG} --cflags ${SBMLPC}`"
        SBML_LIBS=`${PKG_CONFIG} --libs ${SBMLPC}`
      else
        echo "Couldn't find libsbml.pc. Searching for libsbml in candidate directories."
        if test -d /usr/include/sbml
          then SBML_INCDIR="-DHAVE_SBML -I/usr/include/"
          else
            if test -d /usr/local/include/sbml
              then SBML_INCDIR="-DHAVE_SBML -I/usr/local/include/"
            fi;
        fi;
    fi;
    CHECK_XML()
        SBML_INCDIR="$LIBXML_INCDIR $SBML_INCDIR"
        SBML_LIBS="$XML_LIBS"
    CPPFLAGS="${CPPFLAGS} ${SBML_INCDIR}"
    PKG_LIBS="${PKG_LIBS} ${SBML_LIBS} -lsbml"

  else
    echo 'SBML not found.';

  if ${XML}
    then
      echo 'libXML2 found'
      CHECK_XML()
      CPPFLAGS="${CPPFLAGS} ${LIBXML_INCDIR}"
      PKG_LIBS="${PKG_LIBS} ${XML_LIBS}"
    else
      echo 'XML not found.';
  fi;
fi;


dnl subsutitute varialbles in makevars
AC_SUBST(CPPFLAGS)
AC_SUBST(PKG_LIBS)

AC_OUTPUT(src/Makevars)

dnl cleanup after configure
mv src/Makevars.in src/_Makevars.in_;
rm config.*
