#!/usr/bin/env bash

#author : Suzie Tallon 

#date : February 2026

#script to launch the change point analysis recursively on a pair of species for topologies or a species for dS 

# Usage: ./run_recursive.sh dataname
# Exemple : ./run_recursive.sh sin_urm

WLS_SCRIPT="$1"   # ./cp.wls or cp_dS.wls 
SUBSET_DIR="$2"  

run_wls() {
    local dataname="$1"

    echo "=== Running WLS for $dataname ==="
    $WLS_SCRIPT "${SUBSET_DIR}" "$dataname" || true
    echo "=== Finished WLS for $dataname ==="

    local subset_file="${SUBSET_DIR}/${dataname}_subsets.csv"

    # if no subset, no recursion
    if [[ ! -f "${subset_file}" ]]; then
        echo "No subsets for $dataname"
    else
        echo "Subsets found in ${subset_file} — processing recursively..."

        mapfile -t subsets < "${subset_file}"

        for subset_name in "${subsets[@]}"; do
            # ignore empty lines
            [[ -z "$subset_name" ]] && continue

            # remove quote mark if necessary
            subset_name=$(echo "$subset_name" | tr -d '"')

            echo ">>> Processing subset: $subset_name"
            run_wls "$subset_name"
        done
    fi
}

# Point d’entrée
run_wls "$3"
