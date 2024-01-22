import click as ck
import math
import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

class Align(object):

    def __init__(self, chr_id, chr_len, chr_start, chr_end,
                 ctg_id, ctg_len, ctg_start, ctg_end, is_reverse=False):
        self.chr_id = chr_id
        self.chr_len = chr_len
        self.chr_start = chr_start
        self.chr_end = chr_end
        self.ctg_id = ctg_id
        self.ctg_len = ctg_len
        self.ctg_start = ctg_start
        self.ctg_end = ctg_end
        self.is_reverse = is_reverse

    @property
    def ctg_tail(self):
        return self.ctg_len - self.ctg_end

    def toarray(self):
        if self.is_reverse:
            return [f'-{self.ctg_id}', self.ctg_start, self.ctg_end]
        else:
            return [self.ctg_id, self.ctg_start, self.ctg_end]
        
@ck.command()
@ck.option('--alignment-file', '-af', default='verkko_trio_alignment_hap1.txt', help='Filtered alignments from filter_alignment.py')
@ck.option('--output', '-o', default='', help='Output JSON file with configs')
def main(alignment_file, output):
    aligns = {}
    with open(alignment_file) as f:
        for line in f:
            it = line.strip().split()
            ch = it[0]
            ch_len = int(it[1])
            ch_start, ch_end = int(it[2]), int(it[3])
            ctg = it[7]
            ctg_len = int(it[4])
            ctg_start, ctg_end = int(it[5]), int(it[6])
            reverse = it[9] == '-'
            if reverse:
                ctg_start = ctg_len - ctg_end
                ctg_end = ctg_len - int(it[5])
            if ch not in aligns:
                aligns[ch] = []
            aligns[ch].append(Align(
                ch, ch_len, ch_start, ch_end, ctg,
                ctg_len, ctg_start, ctg_end, reverse))

    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')
    chroms.append('chrM')

    res = {}
    for ch in chroms:
        if ch in aligns:
            ctgs = get_contigs(ch, aligns[ch])
            res.update(ctgs)

    with open(output, 'w') as f:
        f.write(json.dumps(res, indent=4))
    
    
def get_contigs(chrom, aligns):
    result = {}
    len_lim = 5000000
    contigs = join_contigs(aligns)
    if contigs[0].ctg_start < len_lim:
        contigs[0].ctg_start = 0
    if contigs[-1].ctg_tail < len_lim:
        contigs[-1].ctg_end = contigs[-1].ctg_len
    for i in range(1, len(contigs)):
        if contigs[i - 1].ctg_id == 'gap':
            if contigs[i].ctg_start < len_lim:
                contigs[i].ctg_start = 0
        if i + 1 < len(contigs) and contigs[i + 1].ctg_id == 'gap':
            if contigs[i].ctg_tail < len_lim:
                contigs[i].ctg_end = contigs[i].ctg_len
    result[chrom] = [ctg.toarray() for ctg in contigs]
    return result


def join_contigs(contigs):
    joined = []
    joined.append(contigs[0])
    ind = 1
    while ind < len(contigs):
        if joined[-1].ctg_id == contigs[ind].ctg_id and joined[-1].is_reverse == contigs[ind].is_reverse:
            if joined[-1].chr_end < contigs[ind].chr_end:
                joined[-1].ctg_end = contigs[ind].ctg_end
                joined[-1].chr_end = contigs[ind].chr_end
        elif joined[-1].chr_end < contigs[ind].chr_end:
            if joined[-1].chr_start <= contigs[ind].chr_start:
                if joined[-1].chr_end > contigs[ind].chr_start:
                    diff = joined[-1].chr_end - contigs[ind].chr_start
                    contigs[ind].chr_start = joined[-1].chr_end
                    contigs[ind].ctg_start += diff
                else:
                    gap_len = 100 # contigs[ind].chr_start - joined[-1].chr_end
                    gap = Align(
                        contigs[ind].chr_id, -1, -1, -1,
                        'gap', gap_len, 0, gap_len, False)
                    joined.append(gap)
                joined.append(contigs[ind])
            else:
                joined[-1] = contigs[ind]
        ind += 1
    return joined



if __name__ == '__main__':
    main()
