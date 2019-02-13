#!/usr/bin/env bash

jupyter nbconvert --ExecutePreprocessor.timeout=-1 --to notebook --execute $1 --output `basename $1`
