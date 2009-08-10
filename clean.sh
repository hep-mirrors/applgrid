#!/bin/sh -f

if [ $BASEDIR != "" ]; then 
  rm -rf $BASEDIR/{lib,include,bld-pdf,bld-nlo}
fi
