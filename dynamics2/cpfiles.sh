
  for yyy in *.c; do

    if  grep "This file is part of yaehmop" $yyy ; then
       cp $yyy yaehmop
    fi

  done
