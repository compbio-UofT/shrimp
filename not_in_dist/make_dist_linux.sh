#!/bin/bash

#get files
pushd `pwd`/../
make clean
popd
sh get_files
git describe --always > ../GIT_VERSION
echo ./GIT_VERSION >> all_files


#find version
version=`cat ../common/version.h  | cut -d "\"" -f 2 | sed 's/\./_/g'`
echo Packaging version $version

#make the 32bit verison 
echo Making 32bit linux...
SHRiMP_LOCAL_FOLDER=SHRiMP_${version}
SHRiMP_FOLDER=`pwd`/${SHRiMP_LOCAL_FOLDER}
TAR_FILENAME=SHRiMP_${version}.lx26.i686.tar
GZ_FILENAME=${TAR_FILENAME}.gz
rm -rf ${SHRiMP_FOLDER}
while read line; do mkdir -p `dirname $(echo ${SHRiMP_FOLDER}/$line | sed 's/\/\.\//\//g')` 2> /dev/null;  ln -s $(dirname `pwd`)/$(echo $line | sed 's/^\.\///g') ${SHRiMP_FOLDER}/$(echo $line | sed 's/^\.\///g'); done < all_files
rm -rf ${SHRiMP_FOLDER}/bin
export CXXFLAGS='-m32 -Kc++ -wd383,981,1572 -axP -O3 -ipo -openmp -DNDEBUG -static-intel'
export CXX=/opt/intel/cc/10.1.015/bin/icc
pushd `pwd`/../
make clean; make
popd
rm -rf ${SHRiMP_FOLDER}/bin
rm -f ${TAR_FILENAME}
tar -hcf ${TAR_FILENAME} ${SHRiMP_LOCAL_FOLDER}
mkdir ${SHRiMP_FOLDER}/bin
cp ../bin/* ${SHRiMP_FOLDER}/bin/
tar -rf ${TAR_FILENAME} ${SHRiMP_LOCAL_FOLDER}/bin
cp ../utils/split-contigs ${SHRiMP_FOLDER}/utils/
tar -rf ${TAR_FILENAME} ${SHRiMP_LOCAL_FOLDER}/utils/split-contigs
rm -f ${GZ_FILENAME}
gzip -9 ${TAR_FILENAME}

#make the 64bit verison 
echo Making 64bit linux...
SHRiMP_LOCAL_FOLDER=SHRiMP_${version}
SHRiMP_FOLDER=`pwd`/${SHRiMP_LOCAL_FOLDER}
TAR_FILENAME=SHRiMP_${version}.lx26.x86_64.tar
GZ_FILENAME=${TAR_FILENAME}.gz
rm -rf ${SHRiMP_FOLDER}
while read line; do mkdir -p `dirname $(echo ${SHRiMP_FOLDER}/$line | sed 's/\/\.\//\//g')` 2> /dev/null;  ln -s $(dirname `pwd`)/$(echo $line | sed 's/^\.\///g') ${SHRiMP_FOLDER}/$(echo $line | sed 's/^\.\///g'); done < all_files
rm -rf ${SHRiMP_FOLDER}/bin
export CXXFLAGS='-m64 -Kc++ -wd383,981,1572 -axP -O3 -ipo -openmp -DNDEBUG -static-intel'
export CXX=/opt/intel/cce/10.1.015/bin/icc
pushd `pwd`/../
make clean; make
popd
rm -rf ${SHRiMP_FOLDER}/bin
rm -f ${TAR_FILENAME}
tar -hcf ${TAR_FILENAME} ${SHRiMP_LOCAL_FOLDER}
mkdir ${SHRiMP_FOLDER}/bin
cp ../bin/* ${SHRiMP_FOLDER}/bin/
tar -rf ${TAR_FILENAME} ${SHRiMP_LOCAL_FOLDER}/bin
cp ../utils/split-contigs ${SHRiMP_FOLDER}/utils/
tar -rf ${TAR_FILENAME} ${SHRiMP_LOCAL_FOLDER}/utils/split-contigs
rm -f ${GZ_FILENAME}
gzip -9 ${TAR_FILENAME}

#make a source version
echo Making source dist...
SHRiMP_LOCAL_FOLDER=SHRiMP_${version}
SHRiMP_FOLDER=`pwd`/${SHRiMP_LOCAL_FOLDER}
TAR_FILENAME=SHRiMP_${version}.src.tar
GZ_FILENAME=${TAR_FILENAME}.gz
rm -rf ${SHRiMP_FOLDER}
while read line; do mkdir -p `dirname $(echo ${SHRiMP_FOLDER}/$line | sed 's/\/\.\//\//g')` 2> /dev/null;  ln -s $(dirname `pwd`)/$(echo $line | sed 's/^\.\///g') ${SHRiMP_FOLDER}/$(echo $line | sed 's/^\.\///g'); done < all_files
rm -rf ${SHRiMP_FOLDER}/bin
pushd `pwd`/../
make clean
popd
rm -rf ${SHRiMP_FOLDER}/bin
mkdir ${SHRiMP_FOLDER}/bin
rm -f ${TAR_FILENAME}
tar -hcf ${TAR_FILENAME} ${SHRiMP_LOCAL_FOLDER}
rm -f ${GZ_FILENAME}
gzip -9 ${TAR_FILENAME}

