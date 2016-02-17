#!/bin/bash
FIXPROJ=fixproj
GENSEED=genseed
PTRACK=PTrack

# Parse configuration file
if [ ! -r $1 ]
then
    echo "No configuration file specified!"
    exit
fi

IFS="="
while read PARAM VALUE
do
    # Strip leading and trailing whitespace, as well as comments
    PARAM=$(expr "$PARAM" : '[ 	]*\([^ 	].*[^ 	]\)')
    VALUE=$(expr "${VALUE%%!*}" : '[ 	]*\([^ 	].*[^ 	]\)')
    
    if [ "$PARAM" == "GRIDFN" ]
    then
        GRIDFN=$VALUE
    fi
    if [ "$PARAM" == "OUTFN" ]
    then
        OUTFN=$VALUE
    fi
    if [ "$PARAM" == "SEEDFN" ]
    then
        SEEDFN=$VALUE
    fi
    if [ "$PARAM" == "PROJ" ]
    then
        PROJ=$VALUE
    fi
    if [ "$PARAM" == "SEEDDAT" ]
    then
        SEEDDAT=$VALUE
    fi
done < $1

# Perform coordinate transformation on the grid file
if [ -n "$PROJ" ]
then
    $FIXPROJ -p "$PROJ" $GRIDFN
    $FIXPROJ -p "$PROJ" -v lonc latc xc yc $GRIDFN
else
    $FIXPROJ -g $GRIDFN
    $FIXPROJ -v lonc latc xc yc $GRIDFN
fi

# If needed, generate the particle seed file
if [ -n "$SEEDDAT" ]
then
    $GENSEED -f $GRIDFN -t -o $SEEDFN $SEEDDAT
fi

# Run the tracking simulation
$PTRACK $1

# Perform coordinate transformation on the simulation output file
$FIXPROJ -i $OUTFN

