import pandas as pd
import numpy as np
import click as ck
from Bio import SeqIO
import os
from pathlib import Path
import gzip

ck.command()
ck.option('--chrom', '-chr', help='Chromosome', default='')
def main(chrom):
    genomes = {
        'chm13': 'reference/chm13v2.0.fa',
        'ksa001.1': 'genomes/ksa001.1v1.0.0.fa',
        'ksa001.2': 'genomes/ksa001.2v1.0.0.fa',
    }
    for i in range(2, 10):
        genomes[f'ksa00{i}.1'] = f'minimap2/ksa00{i}/hifiasm/ksa00{i}.hap1.fa'
        genomes[f'ksa00{i}.2'] = f'minimap2/ksa00{i}/hifiasm/ksa00{i}.hap2.fa'
    if chrom:
        output_file = f'pggb/input_{chrom}.fa'
    else:
        output_file = 'pggb/input.fa'
    with open(output_file, 'w') as w:
        for key, seq_file in genomes.items():
            with open(seq_file, 'r') as f:
                seqs = SeqIO.parse(f, 'fasta')
                for rec in seqs:
                    rec.id = key + '_' + rec.id
                    if chrom:
                        if rec.id == chrom:
                            SeqIO.write(rec, w, 'fasta')
                    else:
                        SeqIO.write(rec, w, 'fasta')
                    
if __name__ == '__main__':
    main()
