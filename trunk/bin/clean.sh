#!/bin/sh -f

if [ "$BASEDIR" = "" ]; then 
  rm -rf lib include bld-pdf bld-nlo
else
  rm -rf $BASEDIR/{lib,include,bld-pdf,bld-nlo}
fi
