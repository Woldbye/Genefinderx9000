#!/bin/bash
IN="${1}"
OUT="${2}"
PROGRAM="${3}"
COUNT=0

printf "<<< Generating an .ans file for each input file with path@'${IN}/*.in' >>>\n"
for F_IN in $IN/*.in 
do
    NAME="$(basename ${F_IN} .in)"
    F_OUT="$OUT/$NAME.ans"
    COUNT=$((COUNT+1))
    printf "[$COUNT] Running test on input@$F_IN"

    # Run program with input and forward output to path@ACT
    $PROGRAM < $F_IN > $F_OUT
    echo " :: SUCCESS => $F_OUT"
done