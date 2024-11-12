import csv
import itertools


def generate_similar_sequences(sequence, replacement_rules_file=None, Replacement_ratio=None):
    """
    Generate similar sequences based on a given sequence and optional replacement rules.

    :param sequence: The input sequence for which similar sequences are to be generated.
    :type sequence: str
    :param replacement_rules_file: (Optional) The file path containing the replacement rules. If not provided, default replacement rules will be used.
    :type replacement_rules_file: str or None
    :param Replacement_ratio: (Optional) The ratio of amino acids to be replaced. If not provided, it will be set to one-third of the total amino acid count.
    :type Replacement_ratio: float or None
    :return: A list of possible similar sequences generated based on the input sequence, replacement rules (either default or from file), and the replacement ratio (either default or specified).
    :rtype: list[str]

    This function generates similar sequences to the given input sequence. 
    """

    # 如果没有指定外部文件，则使用默认规则
    if replacement_rules_file is None:
        replacement_rules = {
            'Glu': ['Asp'],  # 酸性氨基酸
            'Asp': ['Glu'],
            'Lys': ['Arg', 'His'],  # 碱性氨基酸
            'Arg': ['Lys', 'His'],
            'His': ['Lys', 'Arg'],
            'Leu': ['Ile', 'Val', 'Ala'],  # 疏水性氨基酸
            'Ile': ['Leu', 'Val', 'Ala'],
            'Val': ['Leu', 'Ile', 'Ala'],
            'Ala': ['Leu', 'Ile', 'Val'],
            'Phe': ['Tyr', 'Trp'],  # 芳香性氨基酸
            'Tyr': ['Phe', 'Trp'],
            'Trp': ['Phe', 'Tyr'],
            'Ser': ['Thr'],  # 极性氨基酸
            'Thr': ['Ser'],
            'Gly': ['Ala'],  # 小体积氨基酸
            'Cys': ['Ser'],  # 特殊功能基团（极性替换）
            'Met': ['Leu', 'Ile'],  # 疏水/非极性氨基酸
            'Asn': ['Gln'],  # 带酰胺基的极性氨基酸
            'Gln': ['Asn'],
        }
    else:
        replacement_rules = read_replacement_rules_from_file(replacement_rules_file)
    total_aa_count = len(sequence.split('--'))
    # 如果没有指定最大替换数量，则设置为总氨基酸数的三分之一
    if Replacement_ratio is None:
        max_replacements = total_aa_count // 3
    else:
        max_replacements = total_aa_count * Replacement_ratio

    # 使用正则表达式匹配氨基酸和其修饰部分（这里暂不处理修饰部分的特殊逻辑）
    fragments = sequence.split('--')
    all_combinations = []

    # 遍历每一个片段，生成替换选项
    for fragment in fragments:
        aa = fragment
        options = []
        if aa in replacement_rules and replacement_rules[aa]:
            # 随机替换的情况，包括不替换和替换的选项
            options.append(aa)
            for replacement in replacement_rules[aa]:
                options.append(replacement)
        else:
            options.append(aa)

        all_combinations.append(options)

    # 生成所有可能的替换组合
    possible_sequences = []
    for combo in itertools.product(*all_combinations):
        replacements_count = sum(1 for i in range(len(fragments)) if combo[i]!= fragments[i])
        if replacements_count >= 0 and replacements_count <= max_replacements:
            new_sequence = '--'.join(combo)
            possible_sequences.append(new_sequence)

    return possible_sequences


def read_replacement_rules_from_file(file_path):
    """
    Read replacement rules from a specified file.

    :param file_path: The path to the file containing the replacement rules.
    :type file_path: str
    :return: A dictionary where the keys are amino acids and the values are lists of possible replacement amino acids as read from the file.
    :rtype: dict

    This function reads replacement rules from a given file. 
    """

    replacement_rules = {}
    with open(file_path, 'r') as file:
        for line in file.readlines():
            line = line.strip()
            if line:
                key, values_str = line.split(':')
                values = [value.strip() for value in values_str.strip('[]').spanplit(',')]
                replacement_rules[key.strip()] = values
    return replacement_rules


def save_sequences_to_csv(sequences, filename):
    """
    Save a list of sequences to a CSV file.
    
    :param sequences: A list of sequences to be saved to the CSV file.
    :type sequences: list[str]
    :param filename: The name of the CSV file to which the sequences will be saved.
    :type filename: str
    :return: None

    This function takes a list of sequences and saves them to a CSV file with the specified filename. 

    """

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Sequence'])  # 写入标题
        for sequence in sequences:
            writer.writerow([sequence])