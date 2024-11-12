import os
global monomer_path, aas_path, AminoAcids_path, aa_smiles_path

path_of_this_file = os.path.abspath(__file__).replace('/cyclicpeptide/setting.py', '')

monomer_path = os.path.join(path_of_this_file,'states/monomer.tsv')
aas_path = os.path.join(path_of_this_file,'states/aas.txt')
aa_smiles_path = os.path.join(path_of_this_file,'states/aa_smiles.txt')
AminoAcids_path = os.path.join(path_of_this_file,'states/AminoAcids.txt')