#!/bin/sh
# set DATE
if [ -z "$DATE" ]; then
 DATE=$(date "+%x")
fi
DATE=$(date -d $DATE "+%Y-%m-%d")
# set VER
if [ -z "$VER" ]; then
 wget https://cran.r-project.org/package=phreeqc -O index.xhtml 2> /dev/null
 VER=`xmlstarlet sel -t -m "/_:html/_:body/_:table[1]/_:tr[1]/_:td[2]" -v . -n index.xhtml 2> /dev/null`
 ver_major=`echo $VER | cut -d '.' -f 1`
 ver_minor=`echo $VER | cut -d '.' -f 2`
 ver_patch=`echo $VER | cut -d '.' -f 3`
 ver_patch=$((ver_patch+1))
 VER="${ver_major}.${ver_minor}.${ver_patch}"
fi
HEAD=`svn st -v configure.ac | awk '{print $1}')`
if [ "$REL" = 'HEAD' ]; then
 REL="$HEAD"
fi
export VER
# sed files
/bin/sh jenkins-dist.sh -v ${VER} -r ${REL} -d ${DATE} -pr ${TAG}
cd R
make VERSION=${VER} RELEASE_DATE=${DATE} build
make VERSION=${VER} RELEASE_DATE=${DATE} check


