#!/bin/bash

echo 'Building Platypus'
python setup.py build

if [ $? -ne 0 ]
then
    echo ''
    echo 'Setup failed. Check previous lines for errors'
    echo ''
    exit 1
fi

cp -rf build/lib*/* .

if [ $? -ne 0 ]
then
    echo ''
    echo 'Could not copy shared libraries into current directory'
    echo ''
    exit 1
else
    echo ''
    echo ''
    echo 'Finished building Platypus'
    echo ''
    echo ''
fi
