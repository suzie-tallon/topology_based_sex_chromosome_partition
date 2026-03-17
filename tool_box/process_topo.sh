#!/bin/bash

#author : Suzie Tallon 

#date : February 2026

# script to summarize the change point analysis outputs into two single files per pair or species

PAIR="$1"   # ex: lat_dio
DIR="$2"

echo "→ Running for : $PAIR"

FILES=$(ls $DIR/topo*_cutted_"$PAIR"*.csv 2>/dev/null)

if [[ -z "$FILES" ]]; then
    echo "No file topo*_cutted found for pair '$PAIR'"
    exit 1
fi

# outout file
OUT1="$DIR/topo_cutted_${PAIR}_final.csv"
OUT2="$DIR/topo_cutted_${PAIR}_errors_final.csv"

# empty output
> "$OUT1"
> "$OUT2"

VECTOR=""  
FIRST=1

for f in $FILES; do
    echo "  → Reading : $f"

    #second line
    line2=$(sed -n '2p' "$f")

    # add to final vector
    if [[ $FIRST -eq 1 ]]; then
        VECTOR="$line2"
        FIRST=0
    else
        VECTOR="$VECTOR,$line2"
    fi

    # Add lines 3 and above to the error file
    sed -n '3,$p' "$f" >> "$OUT2"
done

# write output
echo "$VECTOR" > "$OUT1"

echo
echo "Done !"
echo "  → $OUT1"
echo "  → $OUT2"

