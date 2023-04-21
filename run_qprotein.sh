#!/bin/bash
# Description: qProtein startup script
# Author: Dou Zhixin

# Setup path
SCRIPT_DIR=$(readlink -f "$(dirname "$0")")
VERSION=$(<$SCRIPT_DIR/version.txt)

DB_BUILDER_PATH="$SCRIPT_DIR/bin/qprotein_db_builder.py"
SEQUENCE_ANNOTATOR_PATH="$SCRIPT_DIR/bin/qprotein_sequence_annotator.py"
STRUCTURE_MAPPER_PATH="$SCRIPT_DIR/bin/qprotein_structure_mapper.py"
conda_activate="$HOME/miniconda3/bin/activate"

# Setup python environment
export PYTHONPATH=$PYTHONPATH:$(pwd)

# ---------------------------------------------------------------------------

# qProtein instructions
usage(){
    cat << EOF

Usage: $0 [command] <options>

makedb
  Required:
    -d, --work_dir <path>           Specify the work directory path
    -D, --db_dir <path>             Specify all of the database directory path
    -f, --fasta_file <path>         Path to the query fasta file

  Optional:
    -q, --only_query <True>         Default False, set True if you need only create the database for query

Common options:
  -v, --version                   Show the version and path infomation
  -h, --help                      Show this help message and exit
EOF
}

# ---------------------------------------------------------------------------

# Define options
OPTIONS=vhqd:D:f:
LONGOPTIONS=version,help,only_query,work_dir:,db_dir:,fasta_file:

ARGS=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS -- "$@")
if [ $? -ne 0 ]; then
    echo "Error: Invalid option or missing argument."
    usage
    exit 1
fi

# Parse options
eval set -- "$ARGS"
while true; do
    case "$1" in
        -v|--version)
            echo "qProtein version: $VERSION"
            echo "qProtein directory: $SCRIPT_DIR"
            exit 0
            ;;
        -h|--help)
            usage
            exit
            ;;
        -q|--only_query)
            only_query=True
            shift
            ;;
        -d|--work_dir)
            work_dir="$2"
            shift 2
            ;;
        -D|--db_dir)
            db_dir="$2"
            shift 2
            ;;
        -f|--fasta_file)
            fasta_file="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "[ERROR] Unknown option: $1"
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------

# Define functions for different modules
function makedb {
    if [ ! -f "$DB_BUILDER_PATH" ]; then
        echo "db_builder script $DB_BUILDER_PATH does not exist."
        exit 1
    fi

    if [[ "$only_query" == "" ]] ; then
        only_query=false
    fi

    if [ "$only_query" == True ] ; then
      if [[ "$work_dir" == "" || "$fasta_file" == "" ]] ; then
        echo "[ERROR] missing required options"
        usage
        exit 1
      fi
    else
      if [[ "$work_dir" == "" || "$db_dir" == "" || "$fasta_file" == "" ]]; then
        echo "[ERROR] missing required options"
        usage
        exit 1
      fi
    fi

    source "$conda_activate" qprotein
    python "$DB_BUILDER_PATH" --work_dir="$1" --db_dir="$2" --fasta_file="$3" --only_query="$4"
}

function sequence_annotator() {
    pass
}

function topt_calculator() {
    pass
}

function structure_mapper() {
    pass
}

function structure_quantizer() {
    pass
}

function results_analyzer() {
    pass
}

# ---------------------------------------------------------------------------

# Select functions
if [ "$1" == "makedb" ]; then
    makedb "$work_dir" "$db_dir" "$fasta_file" "$only_query"
else
    echo "[ERROR] Illegal command received: $1"
    usage
fi

