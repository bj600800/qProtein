from RamachanDraw import plot

pdb_file = r'D:\subject\active\1-PyMulstruct\data\visualization_of_structure_evaluation\wxmy\ranked_0.pdb'

# Drawing the Ramachandran plot
out_file = r'D:\subject\active\1-PyMulstruct\code\tools\structural_analysis\visualization\pamachanFigure.png'
plot(pdb_file, out=out_file)
