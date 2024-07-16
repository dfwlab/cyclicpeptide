import re
from .setting import *

# 打开文本文件并读取内容
file_path = AminoAcids_path  # 替换为你的文件路径
AminoAcids = []
with open(file_path, 'r', encoding='utf-8') as file:
    # 逐行读取文件内容并添加到列表中
    for line in file:
        # 去除行尾的换行符，并将参数添加到列表中
        AminoAcids.append(eval(line.strip()))

################## 序列格式互转 ##################
def read_sequence(sequence):

    """
    **read_sequence** reads sequences in a variety of formats, converting the item into nodes(lists of amino acids) and edges(link information).
    you can also use exact function to transform, such as :py:func:`read_iupac_condensed`, :py:func:`read_graph_representation`, and :py:func:`read_one_letter_sequence`

    :param sequence: cyclic peptide sequence.
    :type sequence: *Format* [Graph presentation, IUPAC condensed, Amino acid chain, One letter code]
    :return: seq_format, nodes, edges
    :rtype: string, List[str], List[(int, int)]

    Example::

        seq_format, nodes, edges = SequenceTransformer.read_sequence('aThr,Tyr,dhAbu,bOH-Gln,Gly,Gln,His,Dab,C13:2(t4.t6)-OH(2.3),Lyx,dhAbu @1,5 @6,10 @0,8')
        print(seq_format, nodes, edges)

    Results::
    
        Graph presentation\n
        ['aThr', 'Tyr', 'dhAbu', 'bOH-Gln', 'Gly', 'Gln', 'His', 'Dab', 'C13:2(t4.t6)-OH(2.3)', 'Lyx', 'dhAbu']\n
        [(1, 5), (6, 10), (0, 8)]

    """

    seq_format = ''
    if '@' in sequence:
        seq_format = 'Graph presentation'
        nodes, edges = read_graph_representation(sequence)
    elif ',' in sequence:
        no_parentheses = re.sub(r'\([^)]*\)', '', sequence)  # 删除所有括号中的逗号
        if ',' in no_parentheses:
            seq_format = 'Graph presentation'
            nodes, edges = read_graph_representation(sequence)
        else:
            seq_format = 'IUPAC condensed'
            nodes, edges = read_iupac_condensed(sequence, sep='-')
    elif '--' in sequence:
        seq_format = 'Amino acid chain'
        nodes, edges = read_iupac_condensed(sequence, sep='--')
    elif '-' in sequence:
        try:
            seq_format = 'IUPAC condensed'
            nodes, edges = read_iupac_condensed(sequence, sep='-')
        except:
            nodes, edges = '', ''
    else:
        seq_format = 'One letter peptide'
        nodes, edges = read_one_letter_sequence(sequence)
    return seq_format, nodes, edges


def create_sequence(nodes, edges):
    """
    **create_sequence** can generate processed nodes(amino acid lists) and edges(link information) into sequences in many different formats.
    you can also use exact function to create, such as :py:func:`create_iupac_condensed`, :py:func:`create_graph_presentation`, :py:func:`create_amino_acid_chain`, and :py:func:`create_one_letter_peptide`

    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: {'iupac_condensed': iupac_condensed, 'amino_acid_chain': amino_acid_chain,
            'graph_presentation': graph_presentation, 'one_letter_peptide': one_letter_peptide}

    Example::

        seq_list = SequenceTransformer.create_sequence(['aThr', 'Tyr', 'dhAbu', 'bOH-Gln', 'Gly', 'Gln', 'His', 'Dab', 'C13:2(t4.t6)-OH(2.3)', 'Lyx', 'dhAbu'], [(1, 5), (6, 10), (0, 8)])
        print(seq_list)
    
    Results::

        {\n
        'iupac_condensed': 'aThr(3)-Tyr(1)-dhAbu-bOH-Gln-Gly-Gln(1)-His(2)-Dab-C13:2(t4.t6)-OH(2.3)(3)-Lyx-dhAbu(2)',\n
        'amino_acid_chain': 'aThr(3)--Tyr(1)--dhAbu--bOH-Gln--Gly--Gln(1)--His(2)--Dab--C13:2(t4.t6)-OH(2.3)(3)--Lyx--dhAbu(2)',\n
        'graph_presentation': 'aThr,Tyr,dhAbu,bOH-Gln,Gly,Gln,His,Dab,C13:2(t4.t6)-OH(2.3),Lyx,dhAbu @1,5 @6,10 @0,8',\n
        'one_letter_peptide': None\n
        }

    """
    try:
        iupac_condensed = create_iupac_condensed(nodes, edges)
    except:
        iupac_condensed = None
    try:
        amino_acid_chain = create_amino_acid_chain(nodes, edges)
    except:
        amino_acid_chain = None
    try:
        graph_presentation = create_graph_presentation(nodes, edges)
    except:
        graph_presentation = None
    try:
        one_letter_peptide = create_one_letter_peptide(nodes)
    except:
        one_letter_peptide = None
    return {'iupac_condensed': iupac_condensed, 'amino_acid_chain': amino_acid_chain,
            'graph_presentation': graph_presentation, 'one_letter_peptide': one_letter_peptide}


def read_iupac_condensed(sequence, sep='-'):
    """

    :param sequence: cyclic peptide sequence.
    :type sequence: IUPAC condensed
    :param sep: Optional "sep".
    :type sep: Amino acid link in sequences
    :return: amino acids, edges
    :rtype: List[str], List[(int, int)]


    """

    # 0. 确认是否有盐或多条链
    if '.' in sequence:
        sequence = sequence.split('.')[0].strip()
    # 1. 确认是否有cyclo[]成环信息
    iscyclo = True if 'cyclo' in sequence else False
    sequence = sequence.replace('cyclo', '').strip('[]')  # 移除 'cyclo' 和方括号
    sequence = re.sub(r'\([^)]*\)', replace_hyphen, sequence)  # 替换（）中的-为～～，（）中为修饰
    items = sequence.split(sep)
    # 2. 确认首尾修饰信息
    header = ''
    tail = ''
    if items[0] in ['H', 'NH2', 'Unk']:
        header = items[0]
        items = items[1:]
    if items[-1] in ['H', 'NH2', 'Unk']:
        tail = items[-1]
        items = items[:-1]
    # 3. 获取氨基酸，DL-修饰添加进氨基酸，首位修饰添加进氨基酸
    amino_acids = []
    edge_marks = []
    pred = ''
    max_edge_mark = 0
    for item in items:
        if item in ['DL', 'D', 'L']:
            pred = item + '-'
            continue
        amino_acid = (pred + item).replace('~~', '-').strip()
        edge_mark = []
        while True:
            edge_mark_match = re.search(r'\(\d+\)$', amino_acid)
            if edge_mark_match:
                amino_acid = amino_acid.replace(edge_mark_match[0], '').strip()
                em = int(edge_mark_match[0][1:-1])
                max_edge_mark = max([max_edge_mark, em])
                edge_mark.append(em)
            else:
                break
        amino_acids.append(amino_acid)
        edge_marks.append(edge_mark)
        pred = ''
    if header:
        amino_acids[0] = header + '-' + amino_acids[0]
    if tail:
        amino_acids[-1] = amino_acids[-1] + '-' + tail
    # 4. 处理最后空格但有edge标签的
    if amino_acids[-1] == '' and edge_marks:
        amino_acids = amino_acids[:-1]
        edge_mark = edge_marks[-1][:]
        edge_marks = edge_marks[:-1]
        edge_marks[-1].extend(edge_mark)
    # 5. 处理氨基酸N(1)标注
    for i, aa in enumerate(amino_acids):
        match = re.search(r'N\((\d)\)', aa)
        if match:
            amino_acids[i] = re.sub(r'N\((\d)\)', 'N-', aa)
            edge_marks[i].append(int(match.group(1)))

    # print(amino_acids)
    # print(edge_marks)
    # 6. edges信息处理
    edges = []
    for i in range(1, max_edge_mark + 1):
        pair = []
        for idx in range(len(edge_marks)):
            if i in edge_marks[idx]:
                pair.append(idx)
        edges.append((pair[0], pair[1]))
    if iscyclo:
        edges.append((0, len(amino_acids) - 1))
    return amino_acids, edges


def read_graph_representation(sequence):
    """

    :param sequence: cyclic peptide sequence.
    :type sequence: graph representation
    :return: amino acids, edges
    :rtype: List[str], List[(int, int)]

    """

    amin_acids = [i.strip() for i in sequence.split(' @')[0].strip().split(',')]
    edges = [tuple([int(j) for j in i.strip().split(',')]) for i in sequence.split('@')[1:]]
    return amin_acids, edges


def read_one_letter_sequence(sequence):
    """

    :param sequence: cyclic peptide sequence.
    :type sequence: one letter code
    :return: amino acids, []
    :rtype: List[str], []

    """
    tail = ''
    if '(NH2)' == sequence.upper()[-5:]:
        sequence = sequence[:-5]
        tail = '(NH2)'
    amino_acid_refs = {i: j for i, j, _, _ in AminoAcids}
    amino_acid_refs['X'] = 'UNK'
    aas = [amino_acid_refs[i] for i in sequence]
    aas[-1] = aas[-1] + tail
    return aas, []


def create_iupac_condensed(nodes, edges):
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: IUPAC condensed
    
    """
    edge_marker = 1
    cyclo = False
    amino_acids = nodes[:]
    for i, j in edges:
        if min(i, j) == 0 and max(i, j) == len(amino_acids) - 1:
            cyclo = True
            continue
        amino_acids[i] += f'({edge_marker})'
        amino_acids[j] += f'({edge_marker})'
        edge_marker += 1
    sequence = '-'.join(amino_acids)
    sequence = f'cyclo[{sequence}]' if cyclo else sequence
    return sequence


def create_amino_acid_chain(nodes, edges):
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: Amino acid chain
    
    """

    edge_marker = 1
    amino_acids = nodes[:]
    for i, j in edges:
        amino_acids[i] += f'({edge_marker})'
        amino_acids[j] += f'({edge_marker})'
        edge_marker += 1
    sequence = '--'.join(amino_acids)
    return sequence


def create_graph_presentation(nodes, edges):
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: Graph presentation
    
    """
    sequence = ','.join(nodes)
    edges_pre = '' if len(edges) == 0 else ' '.join(['@' + ','.join([str(k) for k in pair]) for pair in edges])
    return (sequence + ' ' + edges_pre).strip()


def create_one_letter_peptide(nodes):  # 忽略侧链互作
    """    
    
    :param nodes: lists of amino acids.
    :param edges: link information.
    :type edges: List[(int, int)]
    :return: sequence
    :rtype: one letter peptide
    
    """

    amino_acid_refs = {i: j for j, i, _, _ in AminoAcids}
    amino_acid_refs['UNK'] = 'X'
    tail = ''
    if '(NH2)' == nodes[-1].upper()[-5:]:
        nodes[-1] = nodes[-1][:-5]
        tail = '(NH2)'
    is_essential_aa = [i.upper().strip() in amino_acid_refs.keys() for i in nodes]
    if False in is_essential_aa:
        return None
    else:
        return ''.join([amino_acid_refs[i.upper().strip()] for i in nodes]) + tail


def replace_hyphen(match):
    # 替换括号内的所有 hyphen (-) 为 tilde (~)
    return match.group(0).replace('-', '~~')