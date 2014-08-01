#!/bin/bash

echo 'Building Platypus'
python setup.py build
cp -rf build/lib*/* .

echo ''
echo ''
echo 'Finished building Platypus'
echo ''
echo ''
