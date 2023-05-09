#!/bin/bash
# Description: qProtein startup script
# Author: Dou Zhixin
echo "Current process ID:"  $$

# Setup path for qprotein
SCRIPT_DIR=$(readlink -f "$(dirname "$0")")
VERSION=$(<$SCRIPT_DIR/version.txt)
DB_BUILDER_PATH="$SCRIPT_DIR/bin/qprotein_db_builder.py"
ADD_query_fasta_PATH="$SCRIPT_DIR/bin/qprotein_query_builder.py"
ANNOTATION_PARSER_PATH="$SCRIPT_DIR/bin/qprotein_annotation_parser.py"
STRUCTURE_MAPPER_PATH="$SCRIPT_DIR/bin/qprotein_structure_mapper.py"
STRUCTURE_QUANTIZER_PATH="$SCRIPT_DIR/bin/qprotein_structure_quantizer.py"

# Setup python environment
conda_activate="$HOME/miniconda3/bin/activate"
export PYTHONPATH=$PYTHONPATH:$(pwd)

# ---------------------------------------------------------------------------

# qProtein instructions
function usage(){
    cat << EOF

Usage: $0 [command] <options>
Step 1: makedb    [Warning! CALL ME ONLY at the FIRST TIME!]    Make database for uniprot dat file
    -W, --work_dir <path>           Specify the work directory path for all tasks
    -D, --db_dir <path>             Specify all of the database directory path

Step 2: add_query                           Add query sequences to sql database
    -W, --work_dir <path>           Specify the work directory path for all tasks
    -f, --query_fasta <path>        Path to the query fasta file

Step 3: annotator                           Annotate query sequences with uniprot and functional database
    -W, --work_dir <path>           Specify the work directory path for all tasks
    -D, --db_dir <path>             Specify the root database directory path
    -d, --task_dir <path>           Specific the query directory
    -f, --query_fasta <path>        Path to the query fasta file

Step 4: mapper                              Map protein structures
    -d, --task_dir <path>           Specific the query directory

Step 5: quantizer                           quantize protein structures
    -d, --task_dir <path>           Specific the query directory

Common options:
  -v, --version                     Show the version and path infomation
  -h, --help                        Show this help message and exit
EOF
}

# ---------------------------------------------------------------------------

# Define options
OPTIONS=vhW:D:d:f:
LONGOPTIONS=version,help,work_dir:,db_dir:,task_dir:,query_fasta:

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
        -W|--work_dir)
            work_dir=$(readlink -f "$2")
            shift 2
            ;;
        -D|--db_dir)
            db_dir=$(readlink -f "$2")
            shift 2
            ;;
        -d|--task_dir)
            task_dir=$(readlink -f "$2")
            shift 2
            ;;
        -f|--query_fasta)
            query_fasta="$2"
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
    # setup environment for makedb
    source "$conda_activate" qprotein
    if [ ! -f "$DB_BUILDER_PATH" ]; then
        echo "[ERROR] db_builder script ($DB_BUILDER_PATH) does not exist."
        exit 1
    fi

    # build database for qProtein
    if [[ "$work_dir" == "" || "$db_dir" == "" ]]; then
      echo "[ERROR] missing required options"
      usage
      exit 1
    fi

    python "$DB_BUILDER_PATH" --work_dir="$work_dir" --db_dir="$db_dir"

    # build database for uniprot
    uniprot_fasta=("uniprot_sprot.fasta.gz" "uniprot_trembl.fasta.gz")
    for file in "${uniprot_fasta[@]}";
    do
        file_path="$db_dir/uniprot/sequence/$file"
        echo "$(date +"%Y-%m-%d %H:%M:%S,%3N") [INFO]: Making database using Diamond for $file_path "
        diamond makedb --in $file_path --db "${file_path%%.*}" > /dev/null
        echo "$(date +"%Y-%m-%d %H:%M:%S,%3N") [INFO]: Diamond makedb finished"
    done
#    # build database for merops
#    merops_path="$db_dir/merops/merops.fasta"
#    echo "$(date +"%Y-%m-%d %H:%M:%S,%3N") [INFO]: Making database using Diamond for $merops_path "
#    diamond makedb --in $merops_path --db "$(dirname "$merops_path")/merops"  > /dev/null
#    echo "$(date +"%Y-%m-%d %H:%M:%S,%3N") [INFO]: Diamond makedb finished"
}

function add_query() {
    source "$conda_activate" qprotein
    if [ ! -f "$ADD_query_fasta_PATH" ]; then
        echo "[ERROR] query_builder script ($ADD_query_fasta_PATH) does not exist."
        exit 1
    fi

    if [ ! -f "$query_fasta" ]; then
        echo "[ERROR] query fasta file ($query_fasta) does not exist."
        exit 1
    fi

    if [[ "$work_dir" == "" || "$query_fasta" == "" ]]; then
      echo "[ERROR] missing required options"
      usage
      exit 1
    fi

    python "$ADD_query_fasta_PATH" --work_dir="$work_dir" --query_fasta="$query_fasta"

}

function annotator() {
    source "$conda_activate" qprotein
    chmod u+rwx "$work_dir"

    if [[ "$work_dir" == "" || "$db_dir" == "" || "$task_dir" == "" || "$query_fasta" == "" ]]; then
        echo "[ERROR] missing required options"
        usage
        exit 1
    fi

    if [ ! -f "$ANNOTATION_PARSER_PATH" ]; then
        echo "[ERROR] Sequence annotator script ($ANNOTATION_PARSER_PATH) does not exist."
        exit 1
    fi

    result_dir="$task_dir/annotation_results"

    if [ ! -d "$result_dir" ]; then
      mkdir -p "$result_dir"
      echo "$(date +"%Y-%m-%d %H:%M:%S,%3N") [INFO]: Make directory $result_dir"
    fi

    # uniprot annotation
    uniprot=("uniprot_sprot" "uniprot_trembl")
    for uni in "${uniprot[@]}";
      do
        if [ ! -f "$result_dir/$uni.out" ]; then
          file_path="$db_dir/uniprot/sequence/$uni.dmnd"
          echo "$(date +"%Y-%m-%d %H:%M:%S,%3N") [INFO]: $uni annotation start"
          diamond_output="$result_dir/$uni.out"
          diamond blastp --db "$file_path" --query "$query_fasta" --out "$diamond_output" \
          --outfmt 6 qseqid sseqid pident length qstart qend sstart send qcovhsp evalue \
          --sensitive --max-target-seqs 1 --evalue 1e-5 > /dev/null
          awk -F '\t' '$3>=70 && $4>200 && $9>=90 {print $0}' "$diamond_output" > "${diamond_output%.*}_70_200_90.out"
          echo "$(date +"%Y-%m-%d %H:%M:%S,%3N") [INFO]: $uni annotation finished"
        fi
      done

    # parser annotation
    python $ANNOTATION_PARSER_PATH --work_dir="$work_dir" --task_dir="$task_dir" --query_fasta="$query_fasta"
}

function mapper() {
    source "$conda_activate" qprotein

    if [[ "$task_dir" == "" ]]; then
        echo "[ERROR] missing required options"
        usage
        exit 1
    fi

    if [ ! -f "$STRUCTURE_MAPPER_PATH" ]; then
        echo "[ERROR] db_builder script ($STRUCTURE_MAPPER_PATH) does not exist."
        exit 1
    fi


    python $STRUCTURE_MAPPER_PATH --task_dir "$task_dir"

}

function quantizer() {
    source "$conda_activate" qprotein

    if [[ "$task_dir" == "" ]]; then
        echo "[ERROR] missing required options"
        usage
        exit 1
    fi

    if [ ! -f "$STRUCTURE_QUANTIZER_PATH" ]; then
        echo "[ERROR] db_builder script ($STRUCTURE_QUANTIZER_PATH) does not exist."
        exit 1
    fi

    python $STRUCTURE_QUANTIZER_PATH --task_dir "$task_dir"
}

# ---------------------------------------------------------------------------

# Select functions
if [ "$1" == "makedb" ]; then
    makedb
elif [ "$1" == "annotator" ]; then
    annotator
elif [ "$1" == "add_query" ]; then
    add_query
elif [ "$1" == "mapper" ]; then
    mapper
elif [ "$1" == "quantizer" ]; then
    quantizer
elif [ "$1" == "" ]; then
    echo "[ERROR] Specify a command"
    usage
else
    echo "[ERROR] Got illegal command: $1"
    usage
fi

