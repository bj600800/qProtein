import os.path
import biotite.structure as struc
import biotite.structure.io as strucio



def adjust_residue_numbers(src_file_path, clean_file_path):
    # 读取文件
    structure = strucio.load_structure(src_file_path)
    structure = struc.renumber_res_ids(structure, start=1)
    strucio.save_structure(clean_file_path, structure)


work_dir = r"D:\subject\active\1-qProtein\data\enzymes\GH12\2_StrucMapping\negative"
tar_dir = r"D:\subject\active\1-qProtein\data\enzymes\GH12\2_StrucMapping\test"
file_list = os.listdir(work_dir)
for file in file_list:
    file_path = os.path.join(work_dir, file)
    filename = os.path.split(file)[1]
    clean_file_path = os.path.join(tar_dir, filename)
    adjust_residue_numbers(file_path, clean_file_path)
