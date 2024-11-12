from rdkit import Chem
from rdkit.Chem import AllChem

################## 结构格式互转 ##################
# 自动识别并读人SMILES InChI Molblock SDFblock PDBblock
# https://chemistry.stackexchange.com/questions/34563/pubchem-inchi-smiles-and-uniqueness
# 以InChi=1S/…开头的InChi标识符是标准InChI。在标准InChI中，InChI标识符“对于任何移动氢原子的排列都必须相同”
# 以InChI=1/…开始为非标准InChI，包括一个以/f开头的额外层（fixed-hydrogen layer）。
# 标准InChI生成结构与SMILES不同，而非标准InChI生成结构与SMILES相同，非标准InChI生成生成代码：Chem.MolToInchi(mol1, options='/FixedH')
# 自动识别并读人SMILES InChI Molblock SDFblock PDBblock


# 生成 SMILES InChI InChIKey Molblock SDFblock PDBblock
def output_molecule(mol, pdbblock=None, conformation=None):
    """

    :param mol: mol file.
    :param pdbblock: Optional "conformation".
    :type pdbblock: pdbblock or None
    :param conformation: Optional "conformation".
    :type conformation: conformation or None
    :return: mol

    """

    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchi = Chem.MolToInchi(mol, options='/FixedH')
    inchikey = Chem.InchiToInchiKey(inchi)
    molblock = Chem.MolToMolBlock(mol, includeStereo=True)
    # PDBblock需要通过空间结构保存手性, 因此其他格式转PDB时需要先创建3D构象;
    if pdbblock is None and conformation is not None:
        pdbblock = Chem.MolToPDBBlock(conformation)
    return {'smiles': smiles, 'inchi': inchi, 'inchikey': inchikey, 'molblock': molblock, 'pdbblock': pdbblock}


def predict_3d_conformation(mol):
    """
    Predict the 3D structure of a molecule.

    :param mol: The input molecule for which the 3D conformation is to be predicted.
    :type mol: Chem.Mol
    :return: The molecule with predicted 3D conformation, including added hydrogens and generated 3D coordinates.
    :rtype: Chem.Mol

    This function takes an input molecule and predicts its 3D structure. First, it creates a copy of the input molecule and adds hydrogens to it using `Chem.AddHs`. Then, it generates 3D coordinates for the molecule with added hydrogens by employing `AllChem.EmbedMolecule` with the `AllChem.ETKDG()` method. Finally, the function returns the molecule with the predicted 3D conformation.

    """
    # molecule = Chem.MolFromSmiles(smiles)
    mol_3d = Chem.AddHs(mol.__copy__())  # 添加氢
    AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # 生成 3D 坐标
    return mol_3d


def mol_optimize(mol):
    """
    Optimizes the 3D structure of a molecule using the Universal Force Field (UFF).

    :param mol: The molecule to be optimized. This molecule should have 3D coordinates.
    :type mol: Chem.Mol
    :return: The optimized molecule with updated 3D coordinates.
    :rtype: Chem.Mol

    This function utilizes the Universal Force Field (UFF) to optimize the 3D structure of a given molecule. The input molecule must possess 3D coordinates for the optimization process to be carried out effectively. Once the optimization is complete, the function returns the molecule with its 3D coordinates updated to reflect the optimized structure.

    """

    # 能量最小化：UFF (Universal Force Field)适用于小分子的力场；多肽建模和优化，专门的生物分子建模工具（如 AMBER、CHARMM 或 GROMACS）会更合适
    AllChem.UFFOptimizeMolecule(mol)  # 使用 UFF 力场进行优化
    return mol



def mol2molblock(mol):
    """
    Convert a molecule object to a MolBlock string.

    :param mol: The molecule object to be converted.
    :type mol: Chem.Mol
    :return: The MolBlock string representation of the input molecule.
    :rtype: str

    """

    return Chem.MolToMolBlock(mol)


def mol2pdbblock(mol):
    """
    Convert a molecule object to a PDBBlock string.

    :param mol: The molecule object to be converted.
    :type mol: Chem.Mol
    :return: The PDBBlock string representation of the input molecule.
    :rtype: str

    """

    return Chem.MolToPDBBlock(mol)
