#!/bin/bash
export PYTHONPATH=$PYTHONPATH:$(pwd)

if [ $# -lt 2 ]; then
	echo "Usage: $0 required parameters: work_dir task_name (e.g. bash run_qprotein.sh work_dir sample_name)"
  exit 1
fi

python /home/dzx/qprotein/code/bin/qprotein_structure_processer.py --work_dir "$1" --task_name "$2"
