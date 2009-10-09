#!/bin/sh

# echo "hoppet.sh $1"
# ls -l $1

if [ "$1" == "" ]; then 
  exit -1
fi

if [ "$1" == "--lib" ]; then 
  if [ "$4" == "" ]; then 
      exit -1
  fi
  if [ -e $2 ]; then 
      echo $3 $4  
  fi
fi

if [ -e $1 ]; then 
  echo "-DHOPPET"
fi

