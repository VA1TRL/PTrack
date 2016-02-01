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
$FIXPROJ -p "$PROJ" $GRIDFN
$FIXPROJ -p "$PROJ" -v lonc latc xc yc $GRIDFN

# If needed, generate the particle seed file
if [ -n "$SEEDDAT" ]
then
    $GENSEED -p "$PROJ" -t -o $SEEDFN $SEEDDAT
fi

# Run the tracking simulation
$PTRACK $1

# Perform coordinate transformation on the simulation output file
$FIXPROJ -p "$PROJ" -i $OUTFN

