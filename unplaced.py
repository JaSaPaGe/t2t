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
    chroms.append('chrY')
    chroms.append('chrM')
    
    with open(config) as f:
        json_text = f.read()
        contigs = json.loads(json_text)

    used = set()
    for chrom, vals in contigs.items():
        for contig_id, s, e in vals:
            if contig_id.startswith('-'):
                contig_id = contig_id[1:]
            used.add(contig_id)
    records = {}
    with open(assembly) as f:
        sequences = SeqIO.parse(f, 'fasta')
        for record in sequences:
            if record.id not in used:
                records[record.id] = record
                
    with open(output, 'w') as f:
        for rec_id, rec in records.items():
            SeqIO.write(rec, f, 'fasta')
                
        


if __name__ == '__main__':
    main()
