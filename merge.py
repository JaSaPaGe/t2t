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
    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    records = merge_contigs(assembly, config)

    out_bed = os.path.splitext(output)[0] + '.cen-mask.bed'
    
    with open(output, 'w') as f, open(out_bed, 'w') as cen_mask:
        for chrom in chroms:
            mrec = SeqRecord(Seq(''), id=chrom, name=chrom, description='')
            c_s = 0
            for chunk in ['s', 'c', 'e']:
                rec_id = chrom + '_' + chunk
                if rec_id in records:
                    rec = records[rec_id]
                    if chunk == 's':
                        c_s = len(rec.seq)
                    elif chunk == 'c':
                        c_e = c_s + len(rec.seq)
                    mrec.seq += rec.seq
            cen_mask.write(f'{chrom}\t{c_s}\t{c_e}\n')
            SeqIO.write(mrec, f, 'fasta')
        # SeqIO.write(get_chrM(), f, 'fasta')
        

def get_chrM():
    with open('data/chrMv0.3.0.fa') as f:
        sequences = SeqIO.parse(f, 'fasta')
        rec = next(sequences)
    rec.description = ''
    return rec

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
        
        


if __name__ == '__main__':
    main()
