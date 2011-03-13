#!/bin/sh

if [ -e "hoppet-config" ]; then
   if [ "$1" = "--cxxflags" ]; then
      hopper-config --cxxflags
   elif [ "$1" = "--libs" ]; then
      hoppet-config --libs
   elif [ "$1" = "" ]; then
      echo -DHOPPET
   else
      hoppet-config $1
   fi
fi


