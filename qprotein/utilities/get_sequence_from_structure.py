from Bio import PDB
from Bio.PDB.MMCIFParser import MMCIFParser
import os


def get_seq(struct_path):
	
	parser = PDB.PDBParser()
	structure = parser.get_structure("structure", struct_path)
	model = structure[0]
	sequence = []
	for chain in model:
		for residue in chain:
			if PDB.is_aa(residue):
				sequence.append(PDB.Polypeptide.three_to_one(residue.get_resname()))
			else:
				sequence.append("-")
	
	# Convert the sequence list into a string
	sequence_str = "".join(sequence)
	file_name = os.path.basename(struct_path)
	file_prifx = os.path.splitext(file_name)[0]
	sequence = {file_prifx: sequence_str}
	return sequence


if __name__ == '__main__':
	dir_struct = r"D:\subject\active\2-GetThermo\data\GH5_2\structure"
	output_file = r"D:\subject\active\2-GetThermo\data\GH5_2\gh5_2.fasta"
	with open(output_file, "w") as f:
		file_list = os.listdir(dir_struct)
		for file in file_list:
			print(file)
			file_path = os.path.join(dir_struct, file)
			sequence =get_seq(file_path)
			for key, value in sequence.items():
				f.writelines([">"+key+"\n", value+"\n"])
			
