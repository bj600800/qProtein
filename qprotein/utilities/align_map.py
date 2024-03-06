from Bio import SeqIO


def get_seq_map(align_before, align_after):
	seq_before = [{i.description:i.seq} for i in SeqIO.parse(align_before, "fasta")]
	seq_after = [{i.description:i.seq} for i in SeqIO.parse(align_after, "fasta")]
	seq_dict = {}
	for before in seq_before:
		for k, v in before.items():
			for after in seq_after:
				for ka, va in after.items():
					if k == ka.split(".pdb")[0]:
						seq_dict[k] = (str(v), str(va))
	return seq_dict


def get_id_map(seq_dict):
	id_map = {}
	for name, seq_tuple in seq_dict.items():
		alignment_dict = {}
		original_index = 0
		aligned_index = 0
		original_sequence = seq_tuple[0]
		aligned_sequence = seq_tuple[1]
		while aligned_index < len(aligned_sequence):
			if aligned_sequence[aligned_index] != "-":
				alignment_dict[original_index + 1] = aligned_index + 1
				original_index += 1
		
			aligned_index += 1
		
			if original_index == len(original_sequence):
				break
		id_map[name] = alignment_dict
	return id_map


def run(origi_fasta, align_fasta):
	seq_dict = get_seq_map(origi_fasta, align_fasta)
	id_map = get_id_map(seq_dict)
	return id_map



