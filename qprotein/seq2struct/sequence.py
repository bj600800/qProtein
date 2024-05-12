"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/13

# Description: Process sequences or Uniprot IDs

# ------------------------------------------------------------------------------
"""
import re
from Bio import SeqIO

from qprotein.utilities import logger
logger = logger.setup_log(name=__name__)


def clean_id(identifier):
    special_chars = ["?", "*", "/", ">", "<", "|", "\\", ".", ":"]
    cleaned_id = identifier
    for char in special_chars:
        cleaned_id = cleaned_id.replace(char, "_")
    cleaned_id = re.sub("_+", "_", cleaned_id)
    return cleaned_id

def get_sequence(fasta_file):
    sequences = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        cleaned_id = clean_id(seq_record.id)
        sequences.append([cleaned_id, str(seq_record.seq)])
    # print(sequences)
    logger.info(f"Get sequences: {len(sequences)}")
    return sequences

def get_id(id_file):
    with open(id_file, "r") as f:
        id_list = [line.strip() for line in f]
    # print(id_list)
    logger.info(f"Get uniprot IDs: {len(id_list)}")
    return id_list
    
    
    

