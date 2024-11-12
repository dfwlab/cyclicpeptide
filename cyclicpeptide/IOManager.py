from rdkit import Chem
import os
import pandas as pd
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG, display
import networkx as nx
import matplotlib.pyplot as plt
import io

###############结构输出####################


def save_mol(molblock, output_file):
    """
    Save a molecule represented by a MolBlock string to a file.

    :param molblock: The MolBlock string representation of the molecule to be saved.
    :type molblock: str
    :param output_file: The name of the file to which the molecule will be saved.
    :type output_file: str
    :return: molblock file

    This function takes a MolBlock string (molblock) and saves the corresponding molecule to a file with the specified output_file name. 

    """

    mol = Chem.MolFromMolBlock(molblock)
    if mol is not None:
        Chem.MolToMolFile(mol, output_file)


def save_pdb(pdbblock, output_file):
    """
    Save a molecule represented by a PDBBlock string to a file.

    :param pdbblock: The PDBBlock string representation of the molecule to be saved.
    :type pdbblock: str
    :param output_file: The name of the file to which the molecule will be saved.
    :type output_file: str
    :return: pdbblock file

    This function takes a PDBBlock string (`pdbblock`) and saves the corresponding molecule to a file with the specified `output_file` name.
    """

    mol = Chem.MolFromPDBBlock(pdbblock)
    if mol is not None:
        Chem.MolToPDBFile(mol, output_file)

def save_svg(svg, output_file):
    """
    Save an SVG (Scalable Vector Graphics) string to a file.

    :param svg: The SVG string to be saved.
    :type svg: str
    :param output_file: The name of the file to which the SVG will be saved.
    :type output_file: str
    :return: svg file

    This function takes an SVG string and saves it to a file with the specified output_file name.
    """

    with open(output_file, 'w') as f:
        f.write(svg)

###################图像输出#################
def plot_smiles(smiles, output_file='output.svg', w=600, h=600, isdisplay=False):  # 默认尺寸为600x600
    """
    Plot a molecule represented by its SMILES (Simplified Molecular Input Line Entry System) string and optionally save or display the resulting SVG image.

    :param smiles: The SMILES string representing the molecule to be plotted.
    :type smiles: str
    :param output_file: (Optional) The name of the file to which the SVG image of the molecule will be saved. Defaults to 'output.svg'.
    :type output_file: str
    :param w: (Optional) The width of the SVG image in pixels. Defaults to 600.
    :type w: int
    :param h: (Optional) The height of the SVG image in pixels. Defaults to 600.
    :type h: int
    :param isdisplay: (Optional) A boolean flag indicating whether to display the SVG image. If True, the image will be displayed; if False, only the SVG string will be returned. Defaults to False.
    :type isdisplay: bool
    :return: If `isdisplay` is False, the SVG string representing the plotted molecule is returned. If `isdisplay` is True, nothing is returned as the image is directly displayed.
    :rtype: str or None

    This function first attempts to convert the given SMILES string into a molecule object using `Chem.MolFromSmiles`. 

    """
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        return ''
    # 设置绘图选项，调整分子的大小
    d2d = rdMolDraw2D.MolDraw2DSVG(w, h)
    # 绘制分子
    d2d.DrawMolecule(m)
    d2d.FinishDrawing()
    # 显示图像
    svg = d2d.GetDrawingText()
    if output_file:
        with open(output_file, 'w') as f:
            f.write(svg)
    if isdisplay:
        display(SVG(svg.replace('svg:', '')))
    else:
        return svg


def plot_graph(G, output_file='output.pdf', dpi=60):
    """
    Plot a graph object and optionally save the plot as a PDF file.

    :param G: The graph object to be plotted, typically a networkx graph (e.g., `nx.Graph`).
    :type G: nx.Graph
    :param output_file: (Optional) The name of the file to which the plot will be saved. Defaults to 'output.pdf'.
    :type output_file: str
    :param dpi: (Optional) The dots per inch (resolution) for the saved plot. Defaults to 60.
    :type dpi: int
    :return: pdf file
    
    """

    plt.figure(figsize=(6, 6), dpi=dpi)
    pos = nx.spring_layout(G)  # 生成布局
    nx.draw(G, pos, with_labels=False, node_color='lightblue', node_size=800)
    labels = nx.get_node_attributes(G, 'code')
    # 调整标签位置
    label_pos = {k: [v[0], v[1] + 0] for k, v in pos.items()}
    nx.draw_networkx_labels(G, label_pos, labels=labels)
    if output_file:
        plt.savefig(output_file, format="pdf")
    plt.show()
    

def graph2svg(G, output_file='outpot.svg', dpi=60):
    """
    Convert a graph object to an SVG file and optionally save it.

    :param G: The graph object to be converted, typically a networkx graph (e.g., `nx.Graph`).
    :type G: nx.Graph
    :param output_file: (Optional) The name of the file to which the SVG will be saved. Defaults to 'outpot.svg'.
    :type output_file: str
    :param dpi: (Optional) The dots per inch (resolution) of the graph layout. Defaults to 60.
    :type dpi: int
    :return: svg file

    """

    plt.figure(figsize=(6, 6), dpi=dpi)
    pos = nx.spring_layout(G)  # 生成布局
    svg_buffer = io.StringIO()
    nx.draw(G, pos, with_labels=False, node_color='lightblue', node_size=800)
    labels = nx.get_node_attributes(G, 'code')
    # 调整标签位置
    label_pos = {k: [v[0], v[1] + 0] for k, v in pos.items()}
    nx.draw_networkx_labels(G, label_pos, labels=labels)
    plt.savefig(svg_buffer, format="svg", bbox_inches='tight')
    plt.show()
    if output_file:
        svg = svg_buffer.getvalue()
        with open(output_file, 'w') as f:
            f.write(svg)