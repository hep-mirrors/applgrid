#!/bin/csh

if ( $#argv < 1 ) then 
  echo "need a version argument"
  exit
endif

set VERSION=$1 

cd /tmp

if ( -e applgrid-full-installation-$VERSION ) then 
  echo "directory already exists" 
  exit
endif
   
svn co svn+ssh://svn.hepforge.org/hepforge/svn/applgrid/tags/applgrid-$VERSION applgrid-full-installation-$VERSION 

cd  applgrid-full-installation-$VERSION
rm -rf `find . -name .svn`
cd ..

# full installation
tar -czf applgrid-full-installation-$VERSION.tgz applgrid-full-installation-$VERSION
scp applgrid-full-installation-$VERSION.tgz login.hepforge.org:applgrid/downloads

# applgrid 
cd applgrid-full-installation-$VERSION
cd appl_grid
./configure
make dist

scp applgrid-*.tar.gz login.hepforge.org:applgrid/downloads

