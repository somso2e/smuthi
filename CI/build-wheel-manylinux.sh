#!/bin/bash

# build wheel
for PYBIN in /opt/python/*/bin/; do
#    targetversion = 
  if [[ $PYBIN == *cp27* ]]; then 
    echo "Skipping ${PYBIN}...."
    continue 
  fi
  if [[ $PYBIN == *cp35* ]]; then 
    echo "Skipping ${PYBIN}...."
    continue 
  fi
  echo "Compiling using Python version ${PYBIN}...."
  /opt/usr/bin/python setup.py sdist
	/opt/usr/bin/python setup.py bdist_wheel --bdist-dir=wheelhouse
done

for WHL in wheelhouse/*.whl; do
  auditwheel repair ${WHL} -w dist;
done
