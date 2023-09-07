import click as ck
import math
import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

    def toarray(self):
        if self.is_reverse:
            return [f'-{self.ctg_id}', self.ctg_start, self.ctg_end]
        else:
            return [self.ctg_id, self.ctg_start, self.ctg_end]
        
@ck.command()
@ck.option('--assembly', '-a', default='', help='Assembly contigs')
@ck.option('--config', '-c', default='', help='JSON file with configuration')
@ck.option('--output', '-o', default='', help='Output fasta file')
def main(assembly, config, output):
    # merge(assembly, config, output)
    # merge_hap1()
    contigs_map()
    # alignment()
    # merge_v021()
    # merge_cen()
    # merge()
    # merge_polished()


def merge(assembly, config, output):
    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    records = merge_contigs(assembly, config)
    
    with open(output, 'w') as f:
        for chrom in chroms:
            mrec = SeqRecord(Seq(''), id=chrom, name=chrom)
            for chunk in ['s', 'c', 'e']:
                rec_id = chrom + '_' + chunk
                if rec_id in records:
                    rec = records[rec_id]
                    mrec.seq += rec.seq
            SeqIO.write(mrec, f, 'fasta')
        
    
def merge_contigs(assembly, config):
    with open(config) as f:
        json_text = f.read()
        contigs = json.loads(json_text)

    records = {}
    with open(assembly) as f:
        sequences = SeqIO.parse(f, 'fasta')
        for record in sequences:
            records[record.id] = record

    result = {}
    used = set()
    reversed_contigs = set()
    for ctg, values in contigs.items():
        mrec = SeqRecord(Seq(''), id=ctg, name=ctg)
        for rec_id, start, end in values:
            if rec_id == 'gap':
                mrec.seq += 'N' * 100
                continue
            if rec_id.startswith('-'):
                rec_id = rec_id[1:]
                if rec_id not in reversed_contigs:
                    records[rec_id].seq = records[rec_id].seq.reverse_complement()
                    reversed_contigs.add(rec_id)
            mrec.seq += records[rec_id].seq[start: end]
            used.add(rec_id)
        result[ctg] = mrec
    return result
        

import json

def contigs_map():
    aligns = {}
    with open('verkko_trio_alignment_hap1.txt') as f:
        for line in f:
            it = line.strip().split()
            ch = it[0]
            ch_len = int(it[1])
            ch_start, ch_end = int(it[2]), int(it[3])
            ctg = it[8]
            ctg_len = int(it[4])
            ctg_start, ctg_end = int(it[5]), int(it[6])
            reverse = it[-1] == '-'
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
    chroms.append('chrM')

    res = {}
    for ch in chroms:
        start_id = ch + '_s'
        end_id = ch + '_e'
        if start_id in aligns:
            ctgs = get_contigs(ch, aligns[start_id], aligns[end_id])
        elif end_id in aligns:
            ctgs = get_contigs(ch, None, aligns[end_id])
        res.update(ctgs)
    print(json.dumps(res, indent=4))
    
    
def get_contigs(chrom, start, end):
    cen_contigs = []
    result = {}
    if start is not None:
        start_contigs = join_contigs(start)
        if start_contigs[0].ctg_start < 20000:
            start_contigs[0].ctg_start = 0
        result[chrom + '_s'] =  [ctg.toarray() for ctg in start_contigs]
        ctg = start_contigs[-1]
        if ctg.ctg_end < ctg.ctg_len:
            align = Align(
                chrom + '_c', -1, -1, -1, ctg.ctg_id, ctg.ctg_len, ctg.ctg_end, ctg.ctg_len, ctg.is_reverse)
            cen_contigs.append(align)
    end_contigs = join_contigs(end)
    ctg = end_contigs[0]
    if ctg.ctg_start > 0:
        if len(cen_contigs) > 0 and cen_contigs[0].ctg_id == ctg.ctg_id:
            cen_contigs[0].ctg_end = ctg.ctg_start
        else:
            align = Align(
                chrom + '_c', -1, -1, -1, ctg.ctg_id, ctg.ctg_len, 0, ctg.ctg_start, ctg.is_reverse)
            cen_contigs.append(align)

    result[chrom + '_c'] = [ctg.toarray() for ctg in cen_contigs]
    result[chrom + '_e'] = [ctg.toarray() for ctg in end_contigs]
    return result


def join_contigs(contigs):
    joined = []
    joined.append(contigs[0])
    ind = 1
    while ind < len(contigs):
        if joined[-1].ctg_id == contigs[ind].ctg_id:
            if joined[-1].chr_end < contigs[ind].chr_end:
                joined[-1].ctg_end = contigs[ind].ctg_end
                joined[-1].chr_end = contigs[ind].chr_end
        elif joined[-1].chr_end < contigs[ind].chr_end:
            if joined[-1].chr_start < contigs[ind].chr_start:
                if joined[-1].chr_end > contigs[ind].chr_start:
                    diff = joined[-1].chr_end - contigs[ind].chr_start
                    contigs[ind].chr_start = joined[-1].chr_end
                    contigs[ind].ctg_start += diff
                joined.append(contigs[ind])
            else:
                joined[-1] = contigs[ind]
        ind += 1
    return joined

    
def alignment():
    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrM')

    with open('minimap/verkko_trio_aln_hap2.srt.paf') as f:
        for line in f:
            it = line.strip().split('\t')
            ctg, ctg_len, ctg_start, ctg_end = it[0], int(it[1]), int(it[2]), int(it[3])
            ch, ch_len, ch_start, ch_end = it[5], int(it[6]), int(it[7]), int(it[8])
            match_num = int(it[9])
            base_num = int(it[10])
            identity = round(match_num / base_num * 100, 2)
            seq_len = ch_end - ch_start + 1
            if identity < 90 or seq_len < 10000:
                continue
            # if ch == chrom:
            print(ch, ch_len, ch_start, ch_end, ctg_len, ctg_start, ctg_end, identity, ctg, seq_len, it[4])

    


if __name__ == '__main__':
    main()
