#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File Name: graph_alignment.py
Author: Dingfeng Wu
Creator: Dingfeng Wu
Date Created: 2023-12-27
Last Modified: 2024-01-03
Version: 1.0.0
License: MIT License
Description: The alignment algorithm is developed by CyclicPepedia based on networkx extension, which can convert cyclic peptides into graphical structures and align them accordingly.

Copyright Information: Copyright (c) 2023 dfwlab (https://dfwlab.github.io/)

The code in this script can be used under the MIT License.
"""
import re
from collections import Counter
from itertools import combinations

import networkx as nx
from networkx.algorithms import isomorphism


def sequence_to_node_edge(sequence):
    """
    Example::

        ref_seq = 'Cys(1)(2)--Cys--OH-DL-Val(2)--4OH-Leu--OH-Ile(1)'

        ref_nodes, ref_edges = ga.sequence_to_node_edge(ref_seq)

    :param sequence: amino acid chain format.
    :type sequence: Format conversion by :py:func:`SequenceTransformer.read_sequence`
    :return: Sequence nodes and edges information
    :rtype: Display with :py:func:`create_graph`

    """
    amino_acids = sequence.split('--')
    nodes = []  # 节点信息: [(0, 'Cys'), (1, 'Cys'), ...]
    edges = []  # 边信息: [(0, 1), (1, 2), ...]
    # 用于存储特殊连接信息
    special_connections = {}
    for i, amino_acid in enumerate(amino_acids):
        # 提取氨基酸名称和特殊连接信息
        amino_acid_name = re.sub(r"\(\d+\)", "", amino_acid)
        # print(amino_acid_name)
        if '-' in amino_acid_name:
            items = [k.strip() for k in amino_acid_name.split('-')]
            for j in range(len(items) - 1):
                nodes.append((str(i) + '_' + str(j), items[j].strip().upper()))
                edges.append((str(i) + '_' + str(j), i))
            amino_acid_name = items[-1]
        nodes.append((i, amino_acid_name.strip().capitalize()))
        # 记录普通连接
        if i > 0:
            edges.append((i - 1, i))
        # 处理特殊连接
        special_conn_ids = re.findall(r"\((\d+)\)", amino_acid)
        # print(special_conn_ids)
        for conn_id in special_conn_ids:
            conn_id = int(conn_id) - 1  # 转换为从0开始的索引
            # 存储特殊连接信息
            if conn_id in special_connections:
                edges.append((special_connections[conn_id], i))
            else:
                special_connections[conn_id] = i
    return nodes, edges


def create_graph(nodes, edges):
    """
    Create a graph using the provided nodes and edges.

    :param nodes: A list of tuples, where each tuple contains a node identifier and its corresponding code.
    :type nodes: list[tuple]
    :param edges: A list of tuples representing the edges in the graph. Each tuple contains two node identifiers indicating the connection between nodes.
    :type edges: list[tuple]
    :return: A networkx graph object (`nx.Graph`) constructed with the given nodes and edges.
    :rtype: nx.Graph

    """
    G = nx.Graph()
    for node, code in nodes:
        G.add_node(node, code=code)
    G.add_edges_from(edges)
    return G


def mcs2similarity(ref, mcs):
    """
    Calculate the similarity between a reference graph and the maximum common subgraph (MCS).

    :param ref: The reference graph used for comparison.
    :type ref: nx.Graph
    :param mcs: The maximum common subgraph (MCS) identified between the reference graph and another graph.
    :type mcs: nx.Graph
    :return: The similarity value between the reference graph and the MCS, calculated as the ratio of the total number of nodes and edges in the MCS to the total number of nodes and edges in the reference graph.
    :rtype: float

    This function computes the similarity between a given reference graph (`ref`) and the maximum common subgraph (`mcs`). The similarity is determined by taking the sum of the number of nodes and edges in the MCS and dividing it by the sum of the number of nodes and edges in the reference graph. The resulting value is a float that represents the degree of similarity between the two graphs.
    """

    return float(len(mcs.nodes) + len(mcs.edges)) / (len(ref.nodes) + len(ref.edges))


def mcs_similarity(query, reference):
    """
    Calculate the similarity between a query graph and a reference graph based on the maximum common subgraph (MCS).

    :param query: The query graph for which similarity with the reference graph is to be calculated.
    :type query: nx.Graph
    :param reference: The reference graph to compare against the query graph.
    :type reference: nx.Graph
    :return: A tuple containing two values. The first value is the total number of nodes and edges in the maximum common subgraph (MCS) if it exists; otherwise, it's 0. The second value is the similarity score between the reference graph and the MCS calculated using the `mcs2similarity` function if the MCS exists; otherwise, it's 0.0.
    :rtype: tuple[int, float]

    """

    subgraphs = generate_subgraphs_from_edges(query)
    n, subg, max_subg = find_max_match_subgraph(subgraphs, reference)
    if max_subg:
        return len(max_subg.nodes) + len(max_subg.edges), mcs2similarity(reference, max_subg)
    return 0, 0.0


def graph_similarity(query, reference):
    """

    :param query: amino acid chain format.
    :type query: Format conversion by :py:func:`SequenceTransformer.read_sequence`
    :return: Sequence nodes and edges information
    :rtype: Display with :py:func:`create_graph`

    """
    distance = nx.graph_edit_distance(query, reference, node_match=node_match)
    similarity = 1 - distance / float(len(query.nodes) + len(query.edges) + len(reference.nodes) + len(reference.edges))
    return distance, similarity


def node_match(n1, n2):
    """Node matching function that matches only if the node's code attribute is the same."""
    # print(n1['code'], n2['code'], n1['code'] == n2['code'])
    return n1['code'] == n2['code']


def generate_subgraphs_from_edges(query_G):
    """Generates a subgraph using the given edge set and finds all the connected components."""
    subgraphs = {}
    for n in range(1, len(query_G.edges) + 1):
        for sub_edges in combinations(query_G.edges(), n):
            # Create a subgraph with the selected edges
            subg = query_G.edge_subgraph(sub_edges).copy()
            if nx.is_connected(subg) and subg.edges():
                subgraphs[len(subg.nodes)] = subgraphs.get(len(subg.nodes), [])
                subgraphs[len(subg.nodes)].append(subg)
    return subgraphs


def is_subgraph_of(query, ref):
    """
    Check if a given query graph is a subgraph of a reference graph.

    :param query: The graph that is being checked to see if it is a subgraph of the reference graph.
    :type query: nx.Graph
    :param ref: The reference graph against which the query graph is being compared.
    :type ref: nx.Graph
    :return: True if the query graph is a subgraph of the reference graph (i.e., there is a subgraph in the reference graph that is isomorphic to the query graph), otherwise False.
    :rtype: bool

    """

    GM = isomorphism.GraphMatcher(ref, query, node_match=node_match)
    return GM.subgraph_is_isomorphic()


def find_max_match_subgraph(subgraphs, ref_G):
    """
    Find the maximum match subgraph among a collection of subgraphs with respect to a reference graph.

    :param subgraphs: A dictionary where the keys are some identifiers (presumably related to the subgraphs) and the values are lists of subgraphs. Each subgraph is expected to be of type `nx.Graph`.
    :type subgraphs: dict
    :param ref_G: The reference graph against which the subgraphs are being compared to find the maximum match.
    :type ref_G: nx.Graph
    :return: A tuple containing three elements:
        - The identifier `n` corresponding to the group of subgraphs where the maximum match subgraph was found (if any). If no match is found, this could be an arbitrary value from the `sorted(subgraphs.keys())` iteration.
        - The subgraph `subg` that was being checked when the maximum match subgraph was found (if any). If no match is found, this could be an arbitrary subgraph from the last iteration.
        - The maximum match subgraph `max_subg` if it exists; otherwise, `None`. If a match is found, it will be of type `nx.Graph`.
    :rtype: tuple

    """

    max_subg = None
    for n in sorted(subgraphs.keys()):
        # print(n, len(subgraphs[n]))
        for subg in subgraphs[n]:
            if is_subgraph_of(subg, ref_G):
                max_subg = subg
                break
        else:
            break
    return n, subg, max_subg


def amino_acid_composition_freq(query, ref):
    """
    Calculate the similarity between the amino acid composition frequencies of two sequences (query and reference).

    :param query: The query sequence for which the amino acid composition frequency is to be analyzed.
    :type query: sequence-like (e.g., list, tuple, string)
    :param ref: The reference sequence for comparison.
    :type ref: sequence-like (e.g., list, tuple, string)
    :return: The similarity value between the amino acid composition frequencies of the query and reference sequences. The similarity is calculated as the ratio of the minimum similarity to the maximum similarity. If the maximum similarity is zero, the returned similarity value is zero.
    :rtype: float

    """

    query_freq = Counter(query)
    ref_freq = Counter(ref)
    query_freq_norm = {aa: count / float(len(query_freq)) for aa, count in query_freq.items()}
    ref_freq_norm = {aa: count / float(len(ref_freq)) for aa, count in ref_freq.items()}
    all_amino_acids = set(query_freq_norm.keys()).union(set(ref_freq_norm.keys()))
    # Calculate similarity
    min_similarity = sum(min(query_freq_norm.get(aa, 0), ref_freq_norm.get(aa, 0)) for aa in all_amino_acids)
    max_similarity = sum(max(query_freq_norm.get(aa, 0), ref_freq_norm.get(aa, 0)) for aa in all_amino_acids)
    similarity = min_similarity / max_similarity if max_similarity > 0 else 0
    return similarity
