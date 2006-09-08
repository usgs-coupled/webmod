R=/usr/bin/R
PACKAGE=phreeqc
VERSION=`grep 'Version' ../${PACKAGE}/DESCRIPTION \
         |cut -d':' -f2|tr -d '[:blank:]'`
PACKAGE_DIR=../${PACKAGE}/
PACKAGE_FILE=${PACKAGE}_${VERSION}.tar.gz

default:
	echo "Usage make [build|remove|install|check|clean]"

build:
	${R} CMD build ${PACKAGE_DIR}

remove:
	${R} CMD REMOVE ${PACKAGE}

install: build remove
	${R} CMD INSTALL ${PACKAGE_FILE}

check: build
	${R} CMD check ${PACKAGE_FILE}

clean:
	rm -f ${PACKAGE_FILE}
	rm -rf ${PACKAGE}.Rcheck

