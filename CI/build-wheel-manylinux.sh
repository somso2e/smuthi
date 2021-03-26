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
  ${PYBIN}/python -m pip install numpy==1.12.0
  ${PYBIN}/python setup.py sdist
  ${PYBIN}/python setup.py bdist_wheel
done

mkdir dist_unrepaired
for WHL in dist/*.whl; do
  mv ${WHL} dist_unrepaired
done

for WHL in dist_unrepaired/*.whl; do
  auditwheel repair ${WHL} -w dist
done
