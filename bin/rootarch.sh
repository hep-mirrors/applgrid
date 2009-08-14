#!/bin/sh

if [ -e $ROOTSYS/bin/root-config ]; then 
   m64=`root-config --cflags | grep m64`
   if [ "$m64" != "" ]; then 
     echo "-m64"
   else
     echo "-m32"
   fi
fi 


