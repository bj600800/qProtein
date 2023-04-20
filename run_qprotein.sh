#!/bin/bash
# Description: qProtein startup script
# Author: Dou Zhixin
# ---------------------------------------------------------------------------
# qProtein instructions
usage(){
       echo ""
       echo "Usage: bash $0 <OPTIONS>"
       echo ""
       echo "Required options:"
       echo "-d <work_dir>            Path to directory of working data"
       echo "-n <task_name>           Name of the query task"
       echo ""
       echo "Optional options:"
       echo ""
       exit 1
       }

# Define options
while getopts ":d:n:" i; do
        case "${i}" in
        d)
                work_dir=$OPTARG
        ;;
        n)
                task_name=$OPTARG
        ;;
        esac
done

# Parse input and set defaults
if [[ "$work_dir" == "" || "$task_name" == "" ]] ; then
  usage
fi
# ---------------------------------------------------------------------------
# Setup user's path variables [Change me depend on your installations]
conda_activate="$HOME/miniconda3/bin/activate"

# Setup qProtein config (Temporary)
qprotein_structure_processer_script="$(readlink -f $(dirname "$0"))/bin/qprotein_structure_processer.py"

# Checkout path variables (Temporary)
if [ ! -f "$qprotein_structure_processer_script" ]; then
    echo "qprotein_structure_processer script $qprotein_structure_processer_script does not exist."
    exit 1
fi

# Setup environment
export PYTHONPATH=$PYTHONPATH:$(pwd)
# ---------------------------------------------------------------------------
# Activate conda environment for qProtein
source "$conda_activate" qprotein

# Run qProtein with required parameters
python "$qprotein_structure_processer_script" \
--work_dir="$work_dir" \
--task_name="$task_name"
