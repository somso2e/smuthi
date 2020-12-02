#!/bin/bash

# build wheel
for PYBIN in /opt/python/*/bin/; do
#    targetversion = 
#    if [$PYBIN == **]   continue
    echo "Compiling using pip version ${PYBIN}...."
#    PATH=/opt/usr/bin/:${PYBIN}:${PATH}   pip wheel --no-deps ./ -w wheelhouse
done
