#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File Name: sequence2structure.py
Author: Dingfeng Wu
Creator: Dingfeng Wu
Date Created: 2024-01-03
Last Modified: 2023-01-05
Version: 1.0.0
License: MIT License
Description: Sequence-to-Structure (Seq2struc) is a computing process based on RDkit and the characteristics of cyclic peptide sequences, which can create cyclic peptide sequecne and convert sequence to cyclic peptide SMILES.

Copyright Information: Copyright (c) 2023 dfwlab (https://dfwlab.github.io/)

The code in this script can be used under the MIT License.
"""

import pandas as pd
from IPython.display import SVG, display
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import re
from .setting import *

def seq2stru_essentialAA(sequence, cyclic=True):

    """
    ``"cyclic"`` controls whether it is looped, and when ``"cyclic = True"``, it generates a looped structure, otherwise a chain.
    
    Exmaple1::

        smiles, peptide = sequence2structure.seq2stru_essentialAA(sequence='APG', cyclic=False)
        print('SMILES:', smiles)
        basic.plot_molecule(peptide, w=300, h=300, isdisplay=True)

    SMILES: C[C@H](N)C(=O)N1CCC[C@H]1C(=O)NCC(=O)O
    
    .. image:: s2s.png
        :alt: Description of the image
        :align: center
        :scale: 50%


    :py:func:`basic.plot_molecule` function is used to output the resulting structure, or you can use :py:func:`basic.plot_smiles` .

    Example2::

        smiles, peptide = sequence2structure.seq2stru_essentialAA(sequence='APG', cyclic=True)
        print('SMILES:', smiles)
        basic.plot_smiles(smiles, w=300, h=300, isdisplay=True)
        
    SMILES: C[C@H](N)C(=O)N1CCC[C@H]1C(=O)NCC(=O)O

    .. image:: s2s.png
        :alt: Description of the image
        :align: center
        :scale: 50%

    Enter a sequence in 'Ala-Ala-Cys-Asp' format that will be converted to a single-letter structure by the :py:func:`code2symbol` function.

    Example3::

        smiles, peptide = sequence2structure.seq2stru_essentialAA(sequence='Ala-Ala-Cys-Asp', cyclic=True)
        print('SMILES:', smiles)
        
    SMILES: C[C@@H]1NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CS)NC1=O

    """
    try:
        if '-' in sequence:
            sequence = [code2symbol(code)[0] for code in sequence.split('-')]
            sequence = ''.join(sequence)
        peptide = create_peptide_of_essentialAA(sequence, cyclic=cyclic)
        peptide = Chem.RemoveHs(peptide)  # 移除所有隐式氢原子
        Chem.AssignAtomChiralTagsFromStructure(peptide)  # 重新计算所有隐式氢原子
        smiles = Chem.MolToSmiles(peptide, canonical=True)
        return smiles, peptide
    except:
        return None, None


def seq2stru_no_essentialAA(sequence, references, cyclic=True):
    """
    Example::

        sequence = 'Aad--NMe-Ala--4OH-Thr--3Me-Pro'
        smiles, peptide = sequence2structure.seq2stru_no_essentialAA(sequence, references, cyclic=True)
        print(smiles)
        basic.plot_molecule(peptide, w=300, h=300, isdisplay=True)

    SMILES : CC1CCN2C(=O)C(C(O)CO)NC(=O)C(C)N(C)C(=O)C(CCCC(=O)O)NC(=O)C12

    .. image:: ns2s.png
        :alt: Description of the image
        :align: center
        :scale: 50%

    """
    try:
        sequence = [i.strip() for i in re.sub(r'\(\d+\)', '', sequence).split('--')]
        peptide = create_peptide_of_no_essentialAA(sequence, references, cyclic=cyclic)
        smiles = Chem.MolToSmiles(peptide, canonical=True)
        try:
            peptide = Chem.RemoveHs(peptide)  # 移除所有隐式氢原子
            Chem.AssignAtomChiralTagsFromStructure(peptide)  # 重新计算所有隐式氢原子
            return smiles, peptide
        except:
            return smiles, peptide
    except:
        return None, None


def code2symbol(code):
    """

    :param code: sequence.
    :type code: amino acid chain format
    :return: One letter code

    """

    # 打开文本文件并读取内容
    file_path = AminoAcids_path  # 替换为你的文件路径
    AminoAcids = []
    with open(file_path, 'r', encoding='utf-8') as file:
        """Reads the contents of the file line by line and adds it to the list"""
        for line in file:
            """# 去除行尾的换行符，并将参数添加到列表中"""
            AminoAcids.append(eval(line.strip()))
    c2s = {i[1]: [i[0], i[2]] for i in AminoAcids}
    if code.upper() in c2s.keys():
        return c2s[code.upper()]
    else:
        return None


def link_aa_by_peptide_bond(mol, c_index, n_index):
    """
    Link two amino acids by a peptide bond in a given molecule.

    :param mol: The molecule object representing the amino acids to be linked.
    :type mol: Chem.Mol
    :param c_index: The index of the carbon atom in the molecule where the peptide bond formation will start (usually the carbonyl carbon of one amino acid).
    :type c_index: int
    :param n_index: The index of the nitrogen atom in the molecule where the peptide bond formation will end (usually the amino nitrogen of another amino acid).
    :type n_index: int
    :return: The modified molecule object with the two amino acids linked by a peptide bond.
    :rtype: Chem.Mol

    This function first tries to identify the oxygen and hydrogen atoms of the hydroxyl group attached to the specified carbon atom (`c_index`). It does this by iterating through the neighbors of the carbon atom and looking for an oxygen atom with a single bond to the carbon. Then, it finds the hydrogen atom attached to that oxygen.

    Next, it creates an editable version of the molecule (`emol`) using `Chem.EditableMol`. It then removes the identified hydrogen and oxygen atoms (if they were found) from the editable molecule.

    Finally, it adds a single bond between the specified carbon (`c_index`) and nitrogen (`n_index`) atoms to form the peptide bond. The modified molecule is then retrieved from the editable molecule and returned.
    """

    o_index = None
    h_index = None
    # 查找与该碳原子相连的羟基的氧和氢原子
    for atom in mol.GetAtomWithIdx(c_index).GetNeighbors():
        if atom.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), c_index).GetBondType() == Chem.BondType.SINGLE:
            o_index = atom.GetIdx()
            for h_atom in atom.GetNeighbors():
                if h_atom.GetAtomicNum() == 1:  # 氢原子
                    h_index = h_atom.GetIdx()
    # print(o_index, h_index)
    # 创建一个可编辑的分子
    emol = Chem.EditableMol(mol)
    # 然后再删除羟基氧原子和氢原子
    emol.RemoveAtom(h_index) if h_index else 1
    emol.RemoveAtom(o_index) if o_index else 1
    emol.AddBond(c_index, n_index, order=Chem.rdchem.BondType.SINGLE)
    # 获取修改后的分子
    mol = emol.GetMol()
    return mol

def detect_backbone(mol):
    """
    Detect the backbone of a given molecule.

    :param mol: The molecule object for which the backbone is to be detected.
    :type mol: Chem.Mol
    :return: A tuple containing two elements:

    The first element is the backbone of the molecule as a tuple of atom indices if matches are found; otherwise, an empty string.
    The second element is a list of the backbone atom indices in reverse order if matches are found; otherwise, None.
    
    :rtype: tuple[str, list[int] or None]

    This function aims to identify the backbone of a given molecule. 

    """

    # 识别骨架并将骨架分子编号为0至N，其他侧链分子编号大于N
    # renumber the backbone atoms so the sequence order is correct:
    # N端乙酰化和C端酰胺化: https://zhuanlan.zhihu.com/p/540769025
    # 肽链一般从N端读到C端,C端一般保留一个COOH,因此C端为HO-C(=O),中间一个C,然后加一个肽键NC(=O),后续就都是C和肽键，一直到N端
    # 正反识别都可以，反向就是OC(=O)CNC(=O)CNC(=O)CN...,如果C端酰胺化,则为NC(=O)CNC(=O)CNC(=O)CN...
    # N端不管有没有乙酰化,都会保留N,无乙酰化为-NH2,有乙酰化为-NC(=O)CH3;
    bbsmiles = "C(=O)CN" * len(
        mol.GetSubstructMatches(Chem.MolFromSmiles('NCC=O')))  # 使用GLY最小氨基酸单元识别氨基酸数量，生成骨架generate backbone SMILES
    # backbone = mol.GetSubstructMatches(Chem.MolFromSmiles(bbsmiles))[0] 报错tuple index out of range
    matches = mol.GetSubstructMatches(Chem.MolFromSmiles(bbsmiles))
    backbone = ''
    if matches:  # 检查是否有匹配结果
        backbone = matches[0]
        # 其他操作
    else:
        print("No matches found.")
        return backbone, None # 如果找不到肽键，或者数量不同则退出，比如存在两个线性肽链通过两个链接形成环的情况，需要逐一搜索；
    backbone_idx = list(backbone)
    backbone_idx.reverse()
    return backbone, backbone_idx
 

def create_peptide_of_essentialAA(sequence, cyclic=True):
    """

    :param sequence: sequence of essential AA.
    :type sequence: One letter code
    :param cyclic: Optional "cyclic".
    :type cyclic: True or False
    :return: mol

    """

    try:
        mol = Chem.MolFromSequence(sequence)
    except:
        return None
    if cyclic:
        backbone, backbone_idx = detect_backbone(mol)
        c_index = backbone[0]
        n_index = backbone[-1]
        mol = link_aa_by_peptide_bond(mol, c_index, n_index)
    return mol


def reference_aa_monomer(monomers_path): 
    """

    :param monomers_path: sequence of essential AA.
    :type sequence: One letter code
    :param cyclic: Optional "cyclic".
    :type cyclic: True or False
    :return: mol

    """
    if monomers_path:
        monomers = pd.read_csv(monomers_path, sep='\t', index_col=0)
    else:
        monomers = pd.read_csv(monomer_path, sep='\t', index_col=0)
    references = {}
    for i in monomers.index:
        smile = monomers.loc[i, 'Smiles']
        code = monomers.loc[i, 'Code']
        weight = float(monomers.loc[i, 'Weight'])
        symbol = monomers.loc[i, 'Symbol']
        symbol = code if str(symbol).strip() == '' or str(symbol).upper().strip() == 'NAN' else symbol
        try:
            aa = Chem.MolFromSmiles(smile)
            match = aa.GetSubstructMatches(Chem.MolFromSmiles("NCC=O"))
            if match and match[0]:
                num_atoms = aa.GetNumAtoms()
                references[code] = [aa, symbol, num_atoms, match[0]]
        except:
            pass
    return references

def create_peptide_of_no_essentialAA(sequence, references, cyclic=True):
    """

    :param sequence: sequence of noessential AA.
    :type sequence: One letter code
    :param references: monomers.
    :type references: from :py:func:`reference_aa_monomer`
    :param cyclic: Optional "cyclic".
    :type cyclic: True or False
    :return: mol

    """
    peptide, pep_idx = references[sequence[0]][0], references[sequence[0]][-1]
    pep_c_index, first_n_index = pep_idx[2], pep_idx[0]
    for aa_code in sequence[1:]:
        aa, idx = references[aa_code][0], references[aa_code][-1]
        c_index, n_index = idx[2], idx[0]
        peptide, pep_c_index = connect_two_aa_with_peptide_bond(peptide, pep_c_index, aa, c_index, n_index)
    if cyclic:
        peptide = link_aa_by_peptide_bond(peptide, pep_c_index, first_n_index)
    return peptide

def connect_two_aa_with_peptide_bond(aa1, c_index1, aa2, c_index2, n_index2):
    """
    Connect two amino acids with a peptide bond.

    :param aa1: The first amino acid molecule.
    :type aa1: Chem.Mol
    :param c_index1: The index of the carbon atom in the first amino acid where the peptide bond will be formed.
    :type c_index1: int
    :param aa2: The second amino acid molecule.
    :type aa2: Chem.Mol
    :param c_index2: The index of the carbon atom in the second amino acid which will be involved in the peptide bond formation.
    :type c_index2: int
    :param n_index2: The index of the nitrogen atom in the second amino acid where the peptide bond will be formed.
    :type n_index2: int
    :return: A tuple containing two elements:
        - The connected peptide molecule formed by joining the two amino acids with a peptide bond.
        - The updated index of the carbon atom in the second amino acid (after accounting for the combination of the two amino acids).
    :rtype: tuple[Chem.Mol, int]

    """
    c_index2, n_index2 = c_index2 + aa1.GetNumAtoms() - 1, n_index2 + aa1.GetNumAtoms() - 1
    peptide = Chem.CombineMols(aa1, aa2)
    peptide = link_aa_by_peptide_bond(peptide, c_index1, n_index2)
    return peptide, c_index2


def get_molblock(mol):
    """
    Retrieve the MolBlock representation of a given molecule.

    :param mol: The molecule object for which the MolBlock is to be retrieved.
    :type mol: Chem.Mol
    :return: The MolBlock string representation of the molecule. The `forceV3000=True` parameter indicates that the MolBlock should be in the V3000 format (if available and applicable).
    :rtype: str

    """

    return Chem.MolToMolBlock(mol, forceV3000=True)


def get_molblock_from_smiles(smiles):
    """
    Generate the MolBlock representation of a molecule from its SMILES string.

    :param smiles: The SMILES string representing the molecule.
    :type smiles: str
    :return: The MolBlock string representation of the molecule. The `forceV3000=True` indicates that the MolBlock will be in the V3000 format (if available and applicable).
    :rtype: str

    This function first converts the given SMILES string into a molecule object using `Chem.MolFromSmiles`. If the conversion is successful, it then retrieves the MolBlock representation of the molecule using `Chem.MolToMolBlock` with the `forceV3000=True` option to potentially get the MolBlock in the V3000 format. If the conversion from SMILES to molecule fails, it will likely return `None` as `Chem.MolToMolBlock` will be called on a `None` object (since `mol` would be `None` if the `Chem.MolFromSmiles` conversion fails).
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToMolBlock(mol, forceV3000=True)
