#!/bin/sh

if [ -e "hoppet-config" ]; then
   if [ "$1" = "--cxxflags" ]; then
      hopper-config --cxxflags
   elseif [ "$1" = "--libs" ]; then
      hoppet-config --libs
   elseif [ "$1" = "" ]; then
      echo -DHOPPET
   elseif
      hoppet-config $1
   fi
fi


