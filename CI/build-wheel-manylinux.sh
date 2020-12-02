#!/bin/bash

# build wheel
for PYBIN in /opt/python/*/bin/; do
#    targetversion = 
  if [$PYBIN == *cp27*] then 
    continue 
  fi
  if [$PYBIN == *cp35*] then 
    continue 
  fi
  echo "Compiling using pip version ${PYBIN}...."
#    PATH=/opt/usr/bin/:${PYBIN}:${PATH}   pip wheel --no-deps ./ -w wheelhouse
done
