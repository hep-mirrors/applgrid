#!/bin/csh -f

# set ARCH_TYPE = `/usr/bin/sys`
set ARCH_TYPE = `setarch.sh`

foreach dir ( obj lib exe )
    ./mkdir.sh ./$dir
    ./mkdir.sh ./$dir/$ARCH_TYPE
end




 
