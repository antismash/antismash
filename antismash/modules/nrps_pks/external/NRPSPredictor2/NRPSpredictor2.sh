#!/bin/bash
# adjust base directory
#export NRPS2BASEDIR=


# for JARs
export LIBS=$NRPS2BASEDIR/lib
export BUILD=$NRPS2BASEDIR/build
# for svm models
export DATADIR=$NRPS2BASEDIR/data

java -Ddatadir=$DATADIR -cp $BUILD/NRPSpredictor2.jar:$LIBS/java-getopt-1.0.13.jar:$LIBS/Utilities.jar:$LIBS/libsvm.jar org.roettig.NRPSpredictor2.NRPSpredictor2 $*

