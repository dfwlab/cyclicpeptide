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
    mol = Chem.MolFromMolBlock(molblock)
    if mol is not None:
        Chem.MolToMolFile(mol, output_file)


def save_pdb(pdbblock, output_file):
    mol = Chem.MolFromPDBBlock(pdbblock)
    if mol is not None:
        Chem.MolToPDBFile(mol, output_file)

def save_svg(svg, output_file):
    with open(output_file, 'w') as f:
        f.write(svg)

###################图像输出#################
def plot_smiles(smiles, output_file='output.svg', w=600, h=600, isdisplay=False):  # 默认尺寸为600x600
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


def plot_graph(G, output_file='output.pdf', dpi=60):#ok
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
    

def graph2svg(G, output_file='outpot.svg', dpi=60):#ok
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