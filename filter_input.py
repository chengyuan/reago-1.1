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

import sys
import os

def get_rc(read):
    read = read.upper()
    read = list(read)

    for i in range(len(read)):
        if read[i] not in "ACTGU":
            read[i] = "N"

    read = "".join(read)

    # DNA
    if read.find("T") >= 0:
        alphabet = {'A': 'T',
                    'G': 'C',
                    'C': 'G',
                    'T': 'A',
                    'N': 'N'}
    # RNA
    else:
        alphabet = {'A': 'U',
                    'G': 'C',
                    'C': 'G',
                    'U': 'A',
                    'N': 'N'}
    rev_com = ""
    for char in read:
        rev_com += alphabet[char]

    rev_com = rev_com[::-1]

    return rev_com



def get_fa(fn):
    seq_d = {}
    first = True
    f = open(fn)
    for line in f:
        if line[0] == ">":
            if first != True:
                seq_d[seq_id] = seq.upper()
            else:
                first = False

            seq_id = line.split()[0][1:]
            seq = ""
        else:
            seq += line[:-1]

    seq_d[seq_id] = seq.upper()
    f.close()

    return seq_d, seq_id[-1]


def process_single_file(out_fn):
    d = {}
    start = False
    f = open(out_fn)
    for line in f:
        if line[:15] == "Hit alignments:":
            start = True
        if start == False:
            continue

        if line[:2] == ">>":
            read_id = line.split()[1]
            base = read_id[:-1]
        elif "!" in line:
            data = line.split()
            m_st, m_ed, s_st, s_ed, strand = data[6], data[7], data[9], data[10], data[11]
            if int(s_st) > int(s_ed):
                s_st, s_ed = data[10], data[9]
            d[base] = [m_st, m_ed, s_st, s_ed, strand]
    f.close()
    return d




def process(out_fn_1, out_fn_2):
    d_1 = process_single_file(out_fn_1)
    d_2 = process_single_file(out_fn_2)
    return d_1, d_2



args = sys.argv
try:
    fn_1, fn_2, out_dir, cm_dir, cm, cpu = args[1:]


    if out_dir[:-1] != "/":
        out_dir += "/"

    if cm_dir[:-1] != "/":
        cm_dir += "/"

    if not os.path.exists(fn_1):
        print "Error:", fn_1, "doesn't exist."
        sys.exit(1)


    if not os.path.exists(fn_2):
        print "Error:", fn_2, "doesn't exist."
        sys.exit(1)


    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    print "Indentifying 16S reads"
    if cm == "b":
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "bacteria.cm " + fn_1 + " > " + out_dir + "b_1.out")
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "bacteria.cm " + fn_2 + " > " + out_dir + "b_2.out")
    elif cm == "a":
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "archaea.cm " + fn_1 + " > " + out_dir + "a_1.out")
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "archaea.cm " + fn_2 + " > " + out_dir + "a_2.out")
    elif cm in ["ab", "ba"]:
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "bacteria.cm " + fn_1 + " > " + out_dir + "b_1.out")
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "bacteria.cm " + fn_2 + " > " + out_dir + "b_2.out")
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "archaea.cm " + fn_1 + " > " + out_dir + "a_1.out")
        os.system("cmsearch --cpu " + cpu + " " + cm_dir + "archaea.cm " + fn_2 + " > " + out_dir + "a_2.out")
except:
    print "Usage: python filter_input.py paired_end_1.fasta paired_end_2.fasta output_dir cm_dir cm_to_use num_of_CPU"
    print "paired_end_1.fasta & pared_end_2.fasta: two ends of paired-end file"
    print "output_dir: output directory"
    print "cm_dir: directory containing covariance models for bacteria and archaea"
    print "cm_to_use: b bactera only, a archaea only, ab both"


db_1, end_symbol_1 = get_fa(fn_1)
db_2, end_symbol_2 = get_fa(fn_2)

if "b" in cm:
    b_1, b_2 = process(out_dir + "b_1.out", out_dir + "b_2.out")
else:
    b_1, b_2 = {}, {}

if "a" in cm:
    a_1, a_2 = process(out_dir + "a_1.out", out_dir + "a_2.out")
else:
    a_1, a_2 = {}, {}


b_base = set(b_1.keys()).intersection(set(b_2.keys()))
a_base = set(a_1.keys()).intersection(set(a_2.keys()))



read_cnt = 1
fo = open(out_dir + "filtered.fasta", "w")
included_read = set([])
for base in b_base:
    m_st_1, m_ed_1, s_st_1, s_ed_1, strand_1 = b_1[base]
    m_st_2, m_ed_2, s_st_2, s_ed_2, strand_2 = b_2[base]

    pos_str_1 = " ".join([m_st_1, m_ed_1, s_st_1, s_ed_1])
    pos_str_2 = " ".join([m_st_2, m_ed_2, s_st_2, s_ed_2])

    if strand_1 == strand_2:
        continue

    read_id_1 = base + end_symbol_1
    read_id_2 = base + end_symbol_2

    seq_1 = db_1[read_id_1]
    seq_2 = db_2[read_id_2]

    if strand_1 == "-":
        seq_1 = get_rc(seq_1)
    if strand_2 == "-":
        seq_2 = get_rc(seq_2)

    included_read.add(read_id_1)
    included_read.add(read_id_2)

    fo.write(">" + str(read_cnt) + ".1 " + pos_str_1  + "\n")
    fo.write(seq_1 + "\n")
    fo.write(">" + str(read_cnt) + ".2 "  + pos_str_2  + "\n")
    fo.write(seq_2 + "\n")
    read_cnt += 1


for base in a_base:
    m_st_1, m_ed_1, s_st_1, s_ed_1, strand_1 = a_1[base]
    m_st_2, m_ed_2, s_st_2, s_ed_2, strand_2 = a_2[base]

    pos_str_1 = " ".join([m_st_1, m_ed_1, s_st_1, s_ed_1])
    pos_str_2 = " ".join([m_st_2, m_ed_2, s_st_2, s_ed_2])

    if strand_1 == strand_2:
        continue

    read_id_1 = base + end_symbol_1
    read_id_2 = base + end_symbol_2

    if read_id_1 in included_read:
        continue

    seq_1 = db_1[read_id_1]
    seq_2 = db_2[read_id_2]

    if strand_1 == "-":
        seq_1 = get_rc(seq_1)
    if strand_2 == "-":
        seq_2 = get_rc(seq_2)

    included_read.add(read_id_1)
    included_read.add(read_id_2)

    fo.write(">" + str(read_cnt) + ".1 " + pos_str_1  + "\n")
    fo.write(seq_1 + "\n")
    fo.write(">" + str(read_cnt) + ".2 "  + pos_str_2  + "\n")
    fo.write(seq_2 + "\n")
    read_cnt += 1





fo.close()

