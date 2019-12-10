#!/bin/bash
# Change for prefered install location
idir=/usr/local

# Location of GISS repo
wdir=/usr/local

gts=gts
gerris=gerris
gfsview=gfsview


build_gts=true
build_gerris=true
build_gfsview=true

if $build_gts ; then
    if ( cd $wdir/$gts && sh ./configure --prefix=$idir/local && make -k && make -k install ) \
       > $wdir/build 2>&1 ; then :
    else
	echo
        echo ============ $wdir/$gts: build failed ============
	echo
	cat $wdir/build
	exit 1
    fi
    build_gerris=true
fi


if $build_gerris ; then
    if ( cd $wdir/$gerris && sh ./configure --prefix=$idir/local && make -k && make -k install ) \
       > $wdir/build 2>&1 ; then :
    else
	echo
        echo ============ $wdir/$gerris: build failed ============
	echo
	cat $wdir/build
	exit 1
    fi
    build_gfsview=true
fi


if $build_gfsview ; then
    if ( cd $wdir/$gfsview && sh ./configure --prefix=$idir/local && make -k && make -k install ) \
       > $wdir/build 2>&1 ; then :
    else
	echo
        echo ============ $wdir/$gfsview: build failed ============
	echo
	cat $wdir/build
	exit 1
    fi
fi
