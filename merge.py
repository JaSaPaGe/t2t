import click as ck
import math
import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

@ck.command()
@ck.option('--assembly', '-a', default='', help='Assembly contigs')
@ck.option('--config', '-c', default='', help='JSON file with configuration')
@ck.option('--output', '-o', default='', help='Output fasta file')
def main(assembly, config, output):
    chroms = [f'chr{i}' for i in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrM')
    records = merge_contigs(assembly, config)
    
    with open(output, 'w') as f:
        for chrom in chroms:
            mrec = SeqRecord(Seq(''), id=chrom, name=chrom, description='')
            if chrom in records:
                rec = records[chrom]
                SeqIO.write(rec, f, 'fasta')
        # SeqIO.write(get_chrM(), f, 'fasta')
        

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
                mrec.seq += 'N' * end
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
