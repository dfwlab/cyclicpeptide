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
        :file: Chemical_P.csv
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
    # 理化性质函数
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
          - Indicates non-compliance with Lipinski's rule of five, \n
            suggesting potential issues with bioavailability as an oral drug.
        * - Veber's Rule
          - Shows non-adherence to Veber's rules, \n
            potentially impacting oral bioavailability and permeability.
        * - Ghose Filter
          - A molecular property filter used to assess the \n
            drug-likeness of a compound based on its physicochemical properties.

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
    # 生成RDKit指纹
    fp = RDKFingerprint(mol)
    bit_vector = fp.ToBitString()
    return bit_vector

def cal_daylight_like_fingerprint(mol):
    # 生成拓扑指纹
    fp = rdFingerprintGenerator.GetFPs([mol])
    bit_vector = fp[0].ToBitString()
    return bit_vector


def cal_morgan_fingerprint(mol):
    """ Generate a Morgan fingerprint (radius 2) """
    # 生成Morgan指纹（半径为2）
    morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    bit_vector = morgan_fp.ToBitString()
    return bit_vector


def cal_MACCS_keys(mol):
    # 生成 MACCS keys 指纹
    maccs_fp = MACCSkeys.GenMACCSKeys(mol)
    bit_vector = maccs_fp.ToBitString()
    return bit_vector


def calculate_amino_acid_composition(sequence):
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

