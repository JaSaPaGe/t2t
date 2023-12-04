from Bio import SeqIO
import click as ck
import os
from pathlib import Path


@ck.command()
@ck.option('--genome', '-g', default='', help='Assembled genome')
def main(genome):
    with open(genome) as f:
        sequences = SeqIO.parse(f, 'fasta')

        chroms_dir = Path(os.path.splitext(genome)[0] + '_chroms')
        chroms_dir.mkdir(exist_ok=True)
        for seq in sequences:
            with open(chroms_dir / f'{seq.id}.fa', 'w') as w:
                SeqIO.write(seq, w, 'fasta')

if __name__ == '__main__':
    main()
