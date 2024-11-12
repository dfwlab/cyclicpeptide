from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Lipinski, rdFingerprintGenerator, \
    MACCSkeys, RDKFingerprint
from rdkit.Chem.rdMolDescriptors import CalcTPSA

logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)  # 只显示 CRITICAL 级别的日志

################## 化学性质 ##################
def chemial_physical_properties_from_smiles(smiles):
    """
    Calculate the cyclopeptide chemical property from SMILES.

    :param smiles: Cyclic peptides SMILES information.
    :return: chemial and physical properties
    :rtype: dict{property name: value}

    you can also use :py:func:`cal_chemial_physical_properties` to calculate the properties from mol

    List of the chemical properties:

    .. csv-table:: 
        :file: ../states/Chemical_P.csv
        :header: "Property", "Description"

    """

    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return {}
    properties = cal_chemial_physical_properties(mol)
    properties = cal_rules(properties)
    properties['RDKit_Fingerprint'] = cal_RDKit_fingerprint(mol)
    properties['Daylight_like_Fingerprint'] = cal_daylight_like_fingerprint(mol)
    properties['Morgan_Fingerprint'] = cal_morgan_fingerprint(mol)
    properties['MACCS_Keys'] = cal_MACCS_keys(mol)
    return properties



def cal_chemial_physical_properties(mol):
    """
    Calculate chemical and physical properties of a molecule.

    :param mol: An RDKit molecule object for which properties are to be calculated.
    :type mol: rdkit.Chem.rdchem.Mol
    :return: A dictionary containing the calculated properties. Keys are property names and values are the corresponding calculated values.
    :rtype: dict

    This function computes various chemical and physical properties of a given molecule
    using RDKit's descriptor calculators.
    """

    properties = {
        'Exact_Mass': Descriptors.ExactMolWt(mol),
        'Topological_Polar_Surface_Area': CalcTPSA(mol),
        'Complexity': Descriptors.FpDensityMorgan1(mol),
        'Crippen_LogP': Descriptors.MolLogP(mol),  # 计算LogP Crippen
        'Heavy_Atom_Count': Lipinski.HeavyAtomCount(mol),
        'Hydrogen_Bond_Donor_Count': Lipinski.NumHDonors(mol),
        'Hydrogen_Bond_Acceptor_Count': Lipinski.NumHAcceptors(mol),
        'Rotatable_Bond_Count': Lipinski.NumRotatableBonds(mol),
        'Formal_Charge': Chem.GetFormalCharge(mol),
        'Refractivity': Descriptors.MolMR(mol),
        'Number_of_Rings': rdMolDescriptors.CalcNumRings(mol),
        'Number_of_Atoms': mol.GetNumAtoms(),
    }
    return properties



def cal_rules(properties):
    """
    :param properties: A list of calculated properties.
    :type properties: calculated by :py:func:`chemial_physical_properties_from_smiles`
    :return: properties

    .. list-table::
        :header-rows: 1

        * - Property
          - Description
        * - Rule of Five
          - Indicates non-compliance with Lipinski's rule of five, suggesting potential issues with bioavailability as an oral drug.
        * - Veber's Rule
          - Shows non-adherence to Veber's rules, potentially impacting oral bioavailability and permeability.
        * - Ghose Filter
          - A molecular property filter used to assess the drug-likeness of a compound based on its physicochemical properties.

     """
    
    # Lipinski 规则五
    lipinski_rule_of_five = properties['Hydrogen_Bond_Donor_Count'] <= 5 and \
                            properties['Hydrogen_Bond_Acceptor_Count'] <= 10 and \
                            properties['Exact_Mass'] <= 500 and \
                            properties['Crippen_LogP'] <= 5

    # Veber 规则
    vebers_rule = properties['Topological_Polar_Surface_Area'] <= 140 and \
                  properties['Rotatable_Bond_Count'] <= 10

    # 检查 Ghose Filter 条件
    ghose_filter = 160 <= properties['Exact_Mass'] <= 480 and \
                   0.4 <= properties['Crippen_LogP'] <= 5.6 and \
                   20 <= properties['Number_of_Atoms'] <= 70

    properties['Rule_of_Five'] = lipinski_rule_of_five
    properties["Veber_Rule"] = vebers_rule
    properties['Ghose_Filter'] = ghose_filter
    return properties


################## 分子指纹 ##################

def cal_RDKit_fingerprint(mol):
    """
    Generate an RDKit fingerprint for a given molecule and return it as a bit vector.

    :param mol: The molecule object for which the RDKit fingerprint is to be generated.
    :type mol: Chem.Mol
    :return: The RDKit fingerprint represented as a bit vector (a string of binary digits).
    :rtype: str

    This function takes a molecule object (`mol`) and uses the `RDKFingerprint` function to create an RDKit fingerprint. Then, it converts this fingerprint into a bit string using the `ToBitString` method. The resulting bit vector is then returned, which can be used for various purposes such as similarity calculations or molecule identification in cheminformatics applications.
    """

    # 生成RDKit指纹
    fp = RDKFingerprint(mol)
    bit_vector = fp.ToBitString()
    return bit_vector

def cal_daylight_like_fingerprint(mol):
    """
    Generates a Daylight-type topological fingerprint from a given molecule.
        
    Parameters Mol:
    Mol: RDKit Mol object representing the input molecule.
        
    Return value bit_vector:
    bit_vector: STR type, a Bit string representation of a molecule's topological fingerprint.
    """

    # 生成拓扑指纹
    fp = rdFingerprintGenerator.GetFPs([mol])
    bit_vector = fp[0].ToBitString()
    return bit_vector


def cal_morgan_fingerprint(mol):
    """
    Generate a Morgan fingerprint (radius 2)
    
    This function takes a molecular object as input and generates its Morgan fingerprint.
    The Morgan fingerprint is a type of molecular descriptor used to encode the structural information of a molecule.
    Here, a Morgan fingerprint with a radius of 2 is used, meaning it considers the molecular environment up to 2 bonds away from each atom.
    
    Parameters mol:
    mol - A molecular object, representing the molecule for which to generate the Morgan fingerprint.
    
    Returns bit_vector:
    bit_vector - A string representing the Morgan fingerprint of the molecule.
    """

    # 生成Morgan指纹（半径为2）
    morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    bit_vector = morgan_fp.ToBitString()
    return bit_vector


def cal_MACCS_keys(mol):
    """
    Generate MACCS keys fingerprint for a given molecule.

    This function calculates the MACCS (Molecular ACCess System) keys fingerprint,
    which is a 166-bit structural key descriptor used for molecular similarity searches.

    Parameters:
    mol (rdkit.Chem.rdchem.Mol): An RDKit molecule object for which to generate the MACCS keys fingerprint.

    Returns:
    str: A bit string representation of the MACCS keys fingerprint for the input molecule.
    """
    # 生成 MACCS keys 指纹
    maccs_fp = MACCSkeys.GenMACCSKeys(mol)
    bit_vector = maccs_fp.ToBitString()
    return bit_vector


def calculate_amino_acid_composition(sequence):
    """
    Calculate the composition of amino acids in a given sequence.

    :param sequence: The amino acid sequence for which the composition is to be calculated.
    :type sequence: str
    :return: A list of tuples where each tuple contains an amino acid and its count in the given sequence.
    :rtype: list[tuple[str, int]]

    This function first reads the contents of a text file named 'AminoAcids.txt' (the file path should be adjusted as needed). The file is expected to contain a list of amino acid representations, likely in a specific format. The function reads each line of the file, strips the newline character at the end of each line, and then uses `eval` to convert the string representation back to its original data type (assuming it was originally a tuple or some other data structure). These are then added to the `AminoAcids` list.

    Next, it extracts the first element of each tuple in the `AminoAcids` list to create a list of amino acids (`aa_list`).

    Finally, it iterates through the `aa_list` and for each amino acid, it counts how many times it appears in the given `sequence`. The result is a list of tuples where each tuple contains an amino acid and its count in the sequence. If an amino acid is not present in the sequence, its count will be 0.
    """
    
    # 打开文本文件并读取内容
    file_path = 'AminoAcids.txt'  # 替换为你的文件路径
    AminoAcids = []
    with open(file_path, 'r', encoding='utf-8') as file:
        """Reads the contents of the file line by line and adds it to the list"""
        for line in file:
            """# 去除行尾的换行符，并将参数添加到列表中"""
            AminoAcids.append(eval(line.strip()))

    aa_list = [i[0] for i in AminoAcids]
    counts = [(aa, sequence.count(aa) if aa in sequence else 0) for aa in aa_list]
    return counts


