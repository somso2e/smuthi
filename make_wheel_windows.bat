python setup.py prepare
python setup.py build_ext --inplace --compiler=mingw32 --fcompiler=gnu95 -f
python setup.py bdist_wheel