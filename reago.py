# Software License Agreement (BSD License)
#
# Copyright (c) 2014, Michigan State University
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#  * Neither the name of Willow Garage, Inc. nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import time
import networkx as nx
#import pygraphviz as pgv
import operator


def write_gene(data):
    gene_cnt = 1
    f = open(full_genes_path, "w")
    for path, gene in data:
        f.write(">gene_" + str(gene_cnt) + "_len=" + str(len(gene)) + "\n")
        f.write(gene + "\n")
        gene_cnt += 1
    f.close()


def write_frag(data):
    frag_cnt = 1
    f = open(fragments_path, "w")
    for path, gene in data:
        f.write(">fragment_" + str(frag_cnt) + "_len=" + str(len(gene)) + "\n")
        f.write(gene + "\n")
        frag_cnt += 1
    f.close()


def get_fa(fn):
    r, r_pos, cm_pos = {}, {}, {}
    f = open(fn)
    for line in f:
        if line[0] == ">":
            read_id, m_st, m_ed, s_st, s_ed = line[1:].split()
            r[read_id] = ""
            r_pos[read_id] = [int(s_st), int(s_ed)]
            cm_pos[read_id] = [int(m_st), int(m_ed)]
        else:
            r[read_id] += line[:-1]
    f.close()
    return r, r_pos, cm_pos


def write_fa(sequence_container, filename, width):
    ARROW = ">"
    LINE_BREAK = "\n"

    f = open(filename, "w")
    for elem in sequence_container:
        if type(sequence_container) == list:
            sequence_id, sequence = elem
        elif type(sequence_container) == dict:
            sequence_id, sequence = elem, sequence_container[elem]

        if sequence_id[0] != ARROW:
            sequence_id = ARROW + sequence_id

        f.write(sequence_id + LINE_BREAK)
        if width == 0:
            f.write(sequence + LINE_BREAK)
        else:
            while sequence:
                f.write(sequence[:width] + LINE_BREAK)
                sequence = sequence[width:]
    f.close()


def timestamp():
    return time.asctime()


def n_read_in_node(node):
    read_list = node.split("|")
    return len(read_list)


def initialize_read_pos(read_db):
    read_position = {}
    for read_id in read_db:
        read_position[read_id] = 0
    return read_position


def combine_duplicated_reads(read_db):
    sequence_to_read_id = {}
    for seq_id, seq in read_db.items():
        if seq not in sequence_to_read_id:
            sequence_to_read_id[seq] = []
        sequence_to_read_id[seq].append(seq_id)

    read_db_cleaned = {}
    for seq in sequence_to_read_id:
        new_id = "|".join(sequence_to_read_id[seq])
        read_db_cleaned[new_id] = seq

    return read_db_cleaned


def create_graph_using_rj(read_db, graph_fn):
    G = nx.DiGraph()

    non_dup_fn = rj_dir + graph_fn + ".fasta"
    write_fa(read_db, non_dup_fn, 0)
    os.system("gt readjoiner prefilter -q -des -readset " + rj_dir + graph_fn+ ".set -db " + rj_dir + graph_fn + ".fasta")
    os.system("gt readjoiner overlap -memlimit 100MB -q -l " + str(MIN_OVERLAP) + " -readset " + rj_dir + graph_fn + ".set > /dev/null 2>&1")
    os.system("gt readjoiner spmtest -readset " + rj_dir + graph_fn + ".set.0 -test showlist > " + rj_dir + graph_fn + ".edge.list")

    read_map = {}
    cnt = 0
    f = open(output_dir + "/rj/" + graph_fn + ".set.des")
    for line in f:
        read_map[str(cnt)] = line[:-1]
        cnt += 1
    f.close()

    f = open(output_dir + "/rj/" + graph_fn + ".edge.list")
    for line in f:
        if "-" in line:
            continue
        read_1, read_2, overlap = line.split(" + ")
        read_id_1, read_id_2 = read_map[read_1], read_map[read_2]
        G.add_edge(read_id_1, read_id_2, overlap = int(overlap))
    f.close()

    return G


def correct_sequencing_error(G, ratio):
    if len(G.nodes()) <= 1:
        return

    alignment_to_starting_node = {}
    starting_nodes = []

    visited = set([])
    for node_str in G.nodes():
        if len(G.predecessors(node_str)) == 0:
            starting_nodes.append(node_str) # every node represents a read

    for starting_node in starting_nodes:
        alignment_to_starting_node[starting_node] = {}
        alignment_to_starting_node[starting_node][starting_node], max_st_pos = 0, 0

        queue = [starting_node] # BFS
        while queue != []:
            cur = queue.pop(0)
            if cur in visited:
                continue
            else:
                visited.add(cur) # ownership of cur

            successors = G.successors(cur)
            queue += successors

            for successor in successors:
                overlap = G[cur][successor]['overlap']
                alignment_to_starting_node[starting_node][successor] \
                        = alignment_to_starting_node[starting_node][cur] + READ_LEN - overlap
                max_st_pos = max(max_st_pos, alignment_to_starting_node[starting_node][successor])

        msa = []
        for successor in alignment_to_starting_node[starting_node]:
            msa.append([" " * alignment_to_starting_node[starting_node][successor] + read_db[successor], successor])

        # correcting...
        for i in range(max_st_pos + READ_LEN):
            composition = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0}
            involved_read = []
            for aligned_read, read_id in msa: # ____ACGATGC..ACGATC 23431.CAJ.1
                if i < len(aligned_read) and aligned_read[i] != ' ':
                    composition[aligned_read[i]] += len(read_id.split("|")) + 1
                    involved_read.append(read_id)

            ttl_cnt = sum(composition[k] for k in composition)

            dominant = 'X'
            for base in composition:
                base_cnt = composition[base]
                if float(base_cnt) / ((ttl_cnt - base_cnt) + 1) > ratio:
                    dominant = base
                    break

            if dominant == 'X': # when no dominant base
                continue

            for read_id in involved_read:
                orig_seq = list(read_db[read_id])
                cur_base = orig_seq[i - alignment_to_starting_node[starting_node][read_id]]
                if float(composition[cur_base]) / ttl_cnt < ERROR_CORRECTION_THRESHOLD:
                    orig_seq[i - alignment_to_starting_node[starting_node][read_id]] = dominant
                    read_db[read_id] = "".join(orig_seq)


def correct_sequencing_error_reverse(G, ratio):
    alignment_to_starting_node = {}
    starting_nodes = []

    visited = set([])
    # get starting nodes
    for node_str in G.nodes():
        if len(G.successors(node_str)) == 0:
            starting_nodes.append(node_str) # every node represents a read

    for starting_node in starting_nodes:
        alignment_to_starting_node[starting_node] = {}
        alignment_to_starting_node[starting_node][starting_node], min_st_pos = 0, 0

        queue = [starting_node] # BFS
        while queue != []:
            cur = queue.pop(0)
            if cur in visited:
                continue
            else:
                visited.add(cur) # ownership of cur

            predecessors = G.predecessors(cur)
            queue += predecessors

            for predecessor in predecessors:
                overlap = G[predecessor][cur]['overlap']
                alignment_to_starting_node[starting_node][predecessor] = \
                        alignment_to_starting_node[starting_node][cur] - (READ_LEN - overlap)
                min_st_pos = min(min_st_pos, alignment_to_starting_node[starting_node][predecessor])


        for ancestor in alignment_to_starting_node[starting_node]:
            alignment_to_starting_node[starting_node][ancestor] -= min_st_pos

        align_disp = []
        for ancestor in alignment_to_starting_node[starting_node]:
            align_disp.append([" " * alignment_to_starting_node[starting_node][ancestor] + read_db[ancestor], ancestor])

        # correcting...
        for i in range(-min_st_pos + READ_LEN):
            composition = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            involved_read = []
            for aligned_read, read_id in align_disp: # ____ACGATGC..ACGATC 23431.CAJ.1
                if i < len(aligned_read) and aligned_read[i] != ' ':
                    composition[aligned_read[i]] += len(read_id.split("|")) + 1
                    involved_read.append(read_id)

            ttl_cnt = sum(composition[k] for k in composition)

            dominant = 'X'
            for base in composition:
                base_cnt = composition[base]
                if float(base_cnt) / ((ttl_cnt - base_cnt) + 1) > ratio:
                    dominant = base
                    break

            if dominant == 'X': # when no dominant base
                continue

            for read_id in involved_read:
                orig_seq = list(read_db[read_id])
                cur_base = orig_seq[i - alignment_to_starting_node[starting_node][read_id]]
                if float(composition[cur_base]) / ttl_cnt < ERROR_CORRECTION_THRESHOLD:
                    orig_seq[i - alignment_to_starting_node[starting_node][read_id]] = dominant
                    read_db[read_id] = "".join(orig_seq)


def collapse_graph(G, candidates):
    while True:
        nodes_to_combine = []
        if not candidates:
            all_node = G.nodes()
        else:
            all_node = candidates

        for node in all_node:
            if G.in_degree(node) == 1 and G.out_degree(G.predecessors(node)[0]) == 1:
                nodes_to_combine.append(node)
                if candidates:
                    candidates.remove(node)

        if not nodes_to_combine:
            break

        for node_to_combine in nodes_to_combine:
            predecessor = G.predecessors(node_to_combine)[0]
            predecessors_predecessors = G.predecessors(predecessor)
            successors = G.successors(node_to_combine)

            # update graph
            combined_node = predecessor  + "|" + node_to_combine
            overlap_to_predecessor = G[predecessor][node_to_combine]["overlap"]

            G.add_node(combined_node)
            for predecessors_predecessor in predecessors_predecessors:
                o = G[predecessors_predecessor][predecessor]["overlap"]
                G.add_edge(predecessors_predecessor, combined_node, overlap = o)

            for successor in successors:
                o = G[node_to_combine][successor]["overlap"]
                G.add_edge(combined_node, successor, overlap = o)

            # update sequences
            offset = len(read_db[predecessor]) - overlap_to_predecessor
            for read_id in node_to_combine.split("|"):
                read_position_db[read_id] += offset

            pred_seq = read_db[predecessor]
            node_seq = read_db[node_to_combine]
            combined_seq = pred_seq + node_seq[overlap_to_predecessor:]

            read_db[combined_node] = combined_seq

            # clean up
            G.remove_node(node_to_combine)
            G.remove_node(predecessor)

            del read_db[node_to_combine]
            del read_db[predecessor]

            if node_to_combine in nodes_to_combine:
                nodes_to_combine.remove(node_to_combine)
            if predecessor in nodes_to_combine:
                nodes_to_combine.remove(predecessor)
    return G


def merge_node(src_list, dst, shared, G, direction):
    N_MIS = 3

    dst_seq  = read_db[dst]
    dst_overlap = G[shared][dst]["overlap"]              if direction == 1 else G[dst][shared]["overlap"]
    dst_remaining  = dst_seq[dst_overlap: ]              if direction == 1 else dst_seq[ :-dst_overlap][::-1]

    to_remove = []
    to_merge  = []
    for src in src_list:
        src_seq  = read_db[src]
        src_overlap = G[shared][src]["overlap"]          if direction == 1 else G[src][shared]["overlap"]
        src_remaining  = src_seq[src_overlap: ]          if direction == 1 else src_seq[ :-src_overlap][::-1]

        if n_read_in_node(src) >= 1.2 * n_read_in_node(dst):
            continue

        mis = 0
        for i in range(min(len(src_remaining), len(dst_remaining))):
            if src_remaining[i] != dst_remaining[i]:
                mis += 1
                if mis > N_MIS:
                    break

        if mis > N_MIS:
            if n_read_in_node(src) < TIP_SIZE:
                to_remove.append(src)
            continue

        offset = dst_overlap - src_overlap if direction == 1 else ((len(dst_seq) - dst_overlap) - (len(src_seq) - src_overlap))

        for read_id in src.split("|"):
            read_position_db[read_id] += offset

        to_merge.append(src)

    if not to_remove + to_merge:
        return None

    for n in to_remove:
        G.remove_node(n)

    if to_merge:
        new_node = dst + "|" + "|".join(to_merge)
        G = nx.relabel_nodes(G, {dst: new_node}, copy = False)
        read_db[new_node] = read_db.pop(dst)

        for n in to_merge:
            G.remove_node(n)

        return new_node
    else:
        return dst


def merge_bifurcation(G):
    while True:
        merged = False
        # fork out
        collapse_candidate = set([])
        for node in G.nodes():
            if node not in G.nodes():
                continue

            successors = set(G.successors(node))
            if len(successors) < 2:
                continue

            tip_candidates = set([s for s in successors if G.out_degree(s) == 0])
            if len(tip_candidates) == 0:
                continue

            dst_candidates = successors - tip_candidates
            if len(dst_candidates) == 0:
                dst_node = max([[n_read_in_node(t), t] for t in tip_candidates])[1]
                tip_candidates.remove(dst_node)
            else:
                dst_node = max([[n_read_in_node(d), d] for d in dst_candidates])[1]  # only one dst node
                dst_candidates.remove(dst_node)                                      # remove dst
                dst_candidates = [d for d in dst_candidates if G.out_degree(d) == 0]
                tip_candidates = tip_candidates.union(dst_candidates)

            merged_node = merge_node(tip_candidates, dst_node, node, G, 1)

            if merged_node:
                merged = True
                collapse_candidate.add(node)

        G = collapse_graph(G, list(collapse_candidate))

        # fork in
        collapse_candidate = set([])
        for node in G.nodes():
            if node not in G.nodes():
                continue

            predecessors = set(G.predecessors(node))
            if len(predecessors) < 2:
                continue

            tip_candidates = set([p for p in predecessors if G.in_degree(p) == 0])# and G.out_degree(p) == 1])
            if len(tip_candidates) == 0:
                continue

            dst_candidates = predecessors - tip_candidates
            if len(dst_candidates) == 0:
                dst_node = max([[n_read_in_node(t), t] for t in tip_candidates])[1]
                tip_candidates.remove(dst_node)
            else:
                dst_node = max([[n_read_in_node(d), d] for d in dst_candidates])[1]  # only one dst node
                dst_candidates.remove(dst_node)                                      # remove dst
                dst_candidates = [d for d in dst_candidates if G.in_degree(d) == 0]   # and G.out_degree(d) == 1]  # only if its out-deg is 0, a node will be considered tip
                tip_candidates = tip_candidates.union(dst_candidates)

            merged_node = merge_node(tip_candidates, dst_node, node, G, -1)
            if merged_node:
                merged = True
                collapse_candidate.add(node)

        G = collapse_graph(G, list(collapse_candidate))

        if merged == False:
            break

    G = collapse_graph(G, [])
    return G


def remove_bubble(G):
    while True:
        bubble_removed = False
        all_node = G.nodes()
        collapse_candidate = set([])

        for node in all_node:
            if node not in G.nodes():
                continue

            successors = [s for s in G.successors(node) if G.in_degree(s) == 1 and G.out_degree(s) == 1]
            if len(successors) < 2:
                continue

            d = {}
            for successor in successors:
                to_node = G.successors(successor)[0] # successor has only one successor
                if to_node not in d:
                    d[to_node] = []
                d[to_node].append(successor)

            for to_node in [n for n in d if len(d[n]) > 1]:
                new_node = merge_node(d[to_node][1:], d[to_node][0], node, G, 1)
                if new_node:
                    bubble_removed = True
                    collapse_candidate.add(new_node)

        G = collapse_graph(G, list(collapse_candidate))
        if not bubble_removed:
            break

    return G


def remove_isolated_node(G):
    for node in G.nodes():
        if  not G.in_degree(node) and not G.out_degree(node) and \
            (n_read_in_node(node) < 5 or len(read_db[node]) < READ_LEN * 1.05):
            G.remove_node(node)

    return G


def get_branching_aid(G_orig):
    G = G_orig.reverse(copy = True)
    d = {}
    starting_nodes = []

    for node_str in G.nodes():
        d[node_str] = set([node_str])
        if G.in_degree(node_str) == 0:
            starting_nodes.append(node_str)

    # BFS
    for starting_node in starting_nodes:
        queue = [starting_node]
        while queue != []:
            front = queue.pop(0)
            successors = G.successors(front)
            for successor in successors:
                d[successor] = d[successor].union(d[front])
                if successor not in queue:
                    queue.append(successor)
    return d


def confidence_increment(visited_path, next_node, future_nodes):
    d, weighted_num_pair_end = {}, 0
    for idx, node in enumerate(visited_path):
        for read_id in node.split("|"):
            base, end = read_id.split(".")
            d[base] = len(visited_path) - idx - 1

    for node in future_nodes:
        for read_id in node.split("|"):
            base, end = read_id.split(".")

            if base in d:
                weighted_num_pair_end += 1 * (CONFIDENCE_BASE ** d[base])

    return weighted_num_pair_end


def get_all_path(G, future_nodes, cur_path, paths):
    last_node = cur_path[-1]
    successors = G.successors(last_node)

    # ending node, stop recursion.
    if successors == []:
        paths.append(cur_path)
        return paths
    else:
        if len(successors) > 1:
            candidate = sorted([[confidence_increment(cur_path, s, future_nodes[s]), s] for s in successors])
            next_node = candidate[-1][1]
        else:
            next_node = successors[0]
        return get_all_path(G, future_nodes, cur_path + [next_node], paths)


def get_contig(path, G):
    contig = read_db[path[0]]
    for idx in range(1, len(path)):
        prev, cur = path[idx-1], path[idx]
        seq = read_db[cur]
        overlap = G[prev][cur]["overlap"]
        contig += seq[overlap:]
    return contig


def get_cm_pos(path, contig):
    min_cm_st = float("inf")
    max_cm_ed = 0

    for read_id in [r for r in"|".join(path).split("|") if read_position_db[r] >= 0 and (read_position_db[r] + READ_LEN <= len(contig))]:
        min_cm_st = min(min_cm_st, cm_pos[read_id][0])
        max_cm_ed = max(max_cm_ed, cm_pos[read_id][1])

    return min_cm_st, max_cm_ed


def get_assemblie(G):
    future_nodes = get_branching_aid(G)
    full_genes = []
    scaffold_candidates = []

    starting_nodes = [n for n in G.nodes() if G.in_degree(n) == 0]
    for node in starting_nodes:
        paths = get_all_path(G, future_nodes, [node], [])
        for path in paths:
            contig = get_contig(path, G)
            if len(contig) >= FULL_LENGTH:
                if NEED_DEFLANK:
                    st_pos = min([r_pos[r][0] - read_position_db[r] for r in path[0].split("|")])
                    ed_pos = max([len(read_db_original[r]) - r_pos[r][1] for r in path[-1].split("|")])
                    deflanked_contig = contig[st_pos : len(contig) - ed_pos]
                else:
                    deflanked_contig = contig
                full_genes.append([path, deflanked_contig])
            else:
                m_st, m_ed = get_cm_pos(path, contig)
                if len(contig) > 120:
                    scaffold_candidates.append([path, m_st, m_ed, contig])

    return full_genes, scaffold_candidates


def conf_connect(path_1, path_2):
    d, n_pe = {}, 0
    for idx, node in enumerate(path_1):
        for read_id in node.split("|"):
            base, end = read_id.split(".")
            d[base] = len(path_1) - idx - 1

    for node in path_2:
        for read_id in node.split("|"):
            base, end = read_id.split(".")
            if base in d:
                n_pe += 1 * (CONFIDENCE_BASE ** d[base])
    return n_pe


def calculate_pairwise_segment_confidence(scaffold_candidates):
    n_candidate = len(scaffold_candidates)
    pairwise_confidence = [[0 for i in range(n_candidate)] for j in range(n_candidate)]

    for i, [path_1, m_st_1, m_ed_1, contig_1] in enumerate(scaffold_candidates):
        for j, [path_2, m_st_2, m_ed_2, contig_2] in enumerate(scaffold_candidates):
            if i == j or min(m_ed_1, m_ed_2) - max(m_st_1, m_st_2) < 10:
                pairwise_confidence[i][j] = 0
            else:
                pairwise_confidence[i][j] = max(conf_connect(path_1, path_2), conf_connect(path_2, path_1))

    return pairwise_confidence


def connect_contig(seg_1, m_st_1, m_ed_1, seg_2, m_st_2, m_ed_2):
    if m_st_1 >= m_st_2 and m_ed_1 <= m_ed_2 or m_st_2 >= m_st_1 and m_ed_2 <= m_ed_1:
        return seg_1 if len(seg_1) >= len(seg_2) else seg_2

    if m_st_1 > m_st_2:
        seg_1, seg_2 = seg_2, seg_1

    N_MIS = int(min(len(seg_1), len(seg_2)) * 0.08)
    overlap = 0
    for i in range(min(len(seg_1), len(seg_2)), 10, -1):
        suffix = seg_1[-i:]
        prefix = seg_2[:i]

        n_mis = 0
        for j in range(i):
            if suffix[j] != prefix[j]:
                n_mis += 1
            if n_mis > N_MIS:
                break

        if n_mis <= N_MIS:
            overlap = i
            break

    if overlap > 0:
        return seg_1 + seg_2[overlap:]
    else:
        return seg_1 + "....." + seg_2


def scaffold(scaffold_candidates):
    cont = True
    full_gene = []

    while cont:
        cont = False
        candidate_next = []
        pairwise_confidence = calculate_pairwise_segment_confidence(scaffold_candidates)

        used_candidate = []
        for row_idx, row in enumerate(pairwise_confidence):
            if row_idx in used_candidate:
                continue

            max_conf_idx, max_conf_val = max(enumerate(row), key = operator.itemgetter(1))
            max_conf_idx_rev = max(enumerate(pairwise_confidence[max_conf_idx]), key = operator.itemgetter(1))[0]

            if row_idx == max_conf_idx_rev and max_conf_val != 0:
                used_candidate += [row_idx, max_conf_idx]

                candidate_1 = scaffold_candidates[row_idx]
                candidate_2 = scaffold_candidates[max_conf_idx]
                path_1, m_st_1, m_ed_1, contig_1 = candidate_1
                path_2, m_st_2, m_ed_2, contig_2 = candidate_2

                contig_new = connect_contig(contig_1, m_st_1, m_ed_1, contig_2, m_st_2, m_ed_2)\
                        if m_st_1 < m_st_2 else connect_contig(contig_1, m_st_1, m_ed_1, contig_2, m_st_2, m_ed_2)
                path_new = path_1 + path_2 if m_st_1 < m_st_2 else path_2 + path_1
                m_st_new, m_ed_new = min(m_st_1, m_st_2), max(m_ed_1, m_ed_2)

                if len(contig_new) < FULL_LENGTH:
                    candidate_next.append([path_new, m_st_new, m_ed_new, contig_new])
                else:
                    full_gene.append([path_new, contig_new])

                cont = True
            else:
                candidate_next.append(scaffold_candidates[row_idx])

        if cont:
            scaffold_candidates = candidate_next
        else:
            break

    return full_gene, [[path, contig] for path, du, du, contig in scaffold_candidates]

# for testing purpose
def draw_graph(graph, filename):
    agraph = pgv.AGraph()
    for node in graph.nodes():
        agraph.add_node(node)
    for edge in graph.edges():
        agraph.add_edge(edge)
        node_1, node_2 = edge
        agraph_edge = agraph.get_edge(node_1, node_2)
        agraph_edge.attr["label"] = graph[node_1][node_2]["overlap"]

    agraph.node_attr["shape"] = "box"
    agraph.graph_attr.update(size='80,80')
    agraph.layout()

    agraph.draw(filename, prog = "dot")


def print_help_info():
    print "-----------------------------------------------------"
    print "Usage: python reago.py filename.fasta output_dir -l READ_LENGTH"
    print "Optional parameters:"
    print "-o OVERLAP, default 0.8"
    print "-e ERROR_CORRECTION_THRESHOLD, default 0.05"
    print "-t TIP_SIZE, default 30"
    print "-b PATH_FINDING_PARAMETER, default 10"
    print "Dependencies:"
    print "1. Readjoiner 1.2"
    print "2. Infernal 1.1.1"
    print "-----------------------------------------------------"
    sys.exit(1)




##### setup
arg_range = {"-l" : [1, float("inf"), None], \
                "-o" : [0.5, 1, 0.7], \
                "-e" : [1, float("inf"), 0.05],\
                "-t" : [1, float("inf"), 30],\
                "-b" : [2, 11, 10], \
                "-f" : [1, float("inf"), 1350]}

args = sys.argv
if len(args) < 4:
    print_help_info()

filename = args[1]
output_dir = args[2]
if output_dir[-1] != "/":
    output_dir += "/"

if os.path.exists(filename) == False:
    print "Error: File", "\'" + filename + "\'", "doesn't exist."
    print_help_info()

for idx in range(3, len(args) - 1, 2):
    arg, val = args[idx], float(args[idx+1])
    if arg not in arg_range:
        print "Error - Invalid arg:", arg
        print_help_info()

    min_val, max_val = arg_range[arg][:2]
    if val < min_val or val >= max_val:
        print "Error: Invalid value for", arg
        print "valid range:", "[" + str(min_val) + ", " + str(max_val) + ")"
        print_help_info()
    else:
        arg_range[arg][2] = val

if arg_range["-l"][2] == None:
    print "Error: read length is required"
    print_help_info()

MIN_OVERLAP = int(arg_range["-l"][2] * arg_range["-o"][2])
READ_LEN = int(arg_range["-l"][2])
TIP_SIZE = int(arg_range["-t"][2])
CONFIDENCE_BASE = int(arg_range["-b"][2])
ERROR_CORRECTION_THRESHOLD = float(arg_range["-e"][2])
FULL_LENGTH = int(arg_range["-f"][2])
NEED_DEFLANK = False

graph_path = output_dir + "graph.data"
plot_dir = output_dir + "plot/"
rj_dir = output_dir + "rj/"
full_genes_path = output_dir + "full_genes.fasta"
fragments_path = output_dir + "fragments.fasta"

if os.path.exists(filename) == False:
    print "Error: File", "\'" + filename + "\'", "doesn't exist."
    sys.exit(1)

if os.path.exists(output_dir) == False:
    os.mkdir(output_dir)

if os.path.exists(rj_dir) == False:
    os.mkdir(rj_dir)


##### main procedure starts
print timestamp(), "REAGO (v1.10) started..."
print "Input file:", filename
print "Parameters:"
for arg in arg_range:
    print arg, arg_range[arg][2]

print timestamp(), "Reading input file..."
read_db, r_pos, cm_pos = get_fa(filename)
read_db_original = dict(read_db)
read_position_db = initialize_read_pos(read_db)
read_db = combine_duplicated_reads(read_db)


print timestamp(), "Initializing overlap graph..."
G = create_graph_using_rj(read_db, "graph")
subgraphs = nx.weakly_connected_component_subgraphs(G)

print timestamp(), "Recovering 16S rRNAs..."
full_genes = []
scaffold_candidates = []
for subgraph in subgraphs:
    correct_sequencing_error(subgraph, 5)
    correct_sequencing_error_reverse(subgraph, 5)

    subgraph_read_db = {}
    for node in subgraph.nodes():
        subgraph_read_db[node] = read_db[node]

    subgraph = create_graph_using_rj(subgraph_read_db, "subgraph_temp")
    subgraph = collapse_graph(subgraph, [])
    subgraph = merge_bifurcation(subgraph)
    subgraph = remove_bubble(subgraph)
    subgraph = remove_isolated_node(subgraph)
    subgraph = collapse_graph(subgraph, [])

    full, scaf = get_assemblie(subgraph)
    full_genes += full
    scaffold_candidates += scaf

print timestamp(), "Scaffolding on short 16S rRNA segments..."
scaf, remaining = scaffold(scaffold_candidates)
full_genes += scaf

print timestamp(), "Write to Files..."
write_gene(full_genes)
write_frag(remaining)


print timestamp(), "Done."
print "- Number of 16S rRNAs:", len(full_genes)
print "- Full genes:", full_genes_path
print "- Gene fragments:", fragments_path
