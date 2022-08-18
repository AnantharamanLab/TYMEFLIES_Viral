#! /usr/bin/env python3
# Author: Kristopher Kieft
# University of Wisconsin-Madison
# 2022

import os
import pysam
import subprocess
import argparse
 
class FilterCoverage:
    def __init__(self, sam, bam, out, percent, edit, threads, length, method, interest):
        self.sam = sam
        self.bam = bam
        self.out = out
        self.interest = interest
        self.percent = percent
        self.edit = edit
        self.threads = threads
        self.length = length
        self.method = method
        if not self.out:
            self.set_out()
        
        self.percent = round(1.0 - self.percent,2)
        if self.interest:
            self.get_interest()
        
        if self.sam:
            self.aligned_sam()

        else:
            self.aligned_bam()
        
        self.sort_check()

        if self.method == 'p':
            self.main_percent()
        else:
            self.main_edist()
        
        if self.indexed:
            subprocess.run(f'rm {self.alignment}.bai', shell=True)

    
    def set_out(self):
        if self.method == 'p':
            base = f'pid{int(self.percent*100)}'
        else:
            base = f'dist{self.edit}'
        if self.sam:
            self.out = self.sam.rsplit('.',1)[0] + f'.{base}.sorted.bam'
        elif self.bam:
            self.out = self.bam.rsplit('.',1)[0] + f'.{base}.sorted.bam'
        if os.path.exists(self.out):
            print('\nError: output filtered BAM already exists. Exiting.')
            print(f'{self.out}\n')
            exit()
    
    def get_interest(self):
        with open(self.interest) as f:
            self.names = set(f.read().split('\n'))
    
    def aligned_sam(self):
        check_aligned = subprocess.check_output(f"samtools view -@ {self.threads} {self.sam} | head -n 1", shell=True)
        if len(check_aligned) == 0:
            print(f'\nNo aligned reads in {self.sam}.\n')
            exit()
        try:
            temp = self.sam.rsplit('/',1)[1]
            self.bam = temp.rsplit('.',1)[0] + '.bam'
        except IndexError:
            self.bam = self.sam.rsplit('.',1)[0] + '.bam'
        if os.path.exists(self.bam):
            print('\nError: looks like a BAM file already exists. Exiting.')
            print(f'{self.bam}\n')
            exit()
        subprocess.run(f"samtools view -@ {self.threads} -S -b {self.sam} > {self.bam}", shell=True)
    
    def aligned_bam(self):
        check_aligned = subprocess.check_output(f"samtools view -@ {self.threads} {self.bam} | head -n 1", shell=True)
        if len(check_aligned) == 0:
            print(f'\nNo aligned reads in {self.bam}.\n')
            exit()

    def sort_check(self):
        sort_check = False
        try:
            check = subprocess.check_output(f"samtools view -@ {self.threads} -H {self.bam} | grep '@HD'", shell=True)
            if "coordinate" in str(check):
                sort_check = True
        except Exception:
            # no @HD line
            pass
        if not sort_check:
            sbam = self.bam.rsplit('.',1)[0] + '.sorted.bam'
            if os.path.exists(sbam):
                sbam = self.bam.rsplit('.',1)[0] + '.py-sorted.bam'

            subprocess.run(f"samtools sort -@ {self.threads} -o {sbam} {self.bam}", shell=True)
            self.alignment = sbam
        else:
            self.alignment = self.bam

        self.indexed = False
        if not os.path.exists(f'{self.alignment}.bai'):
            self.indexed = True
            subprocess.run(f'samtools index {self.alignment}', shell=True)

    def align_id(self, ed, rl, x):
        if ed/rl <= self.percent and rl >= self.length:
            self.outfile.write(x)

    def align_e(self, ed, rl, x):
        if ed <= self.edit and rl >= self.length:
            self.outfile.write(x)

    def main_percent(self):

        bamfile = pysam.AlignmentFile(self.alignment, "rb")
        self.outfile = pysam.AlignmentFile(self.out, "wb", template=bamfile)

        if self.interest:
            for x in bamfile.fetch(until_eof=True):
                genome = x.reference_name
                if genome in self.names:
                    rl = x.query_length
                    ed = 0
                    for t in x.tags:
                        if t[0] == 'NM':
                            ed = t[1]
                            self.align_id(ed, rl, x)
                            break

        else:
            for x in bamfile.fetch(until_eof=True):
                rl = x.query_length
                ed = 0
                for t in x.tags:
                    if t[0] == 'NM':
                        ed = t[1]
                        self.align_id(ed, rl, x)
                        break
    
        bamfile.close()
        self.outfile.close()
    
    def main_edist(self):

        bamfile = pysam.AlignmentFile(self.alignment, "rb")
        self.outfile = pysam.AlignmentFile(self.out, "wb", template=self.alignment)

        if self.interest:
            for x in bamfile.fetch(until_eof=True):
                genome = x.reference_name
                if genome in self.names:
                    rl = x.query_length
                    ed = False
                    for t in x.tags:
                        if t[0] == 'NM':
                            ed = t[1]
                            self.align_e(ed, rl, x)
                            break

        else:
            for x in bamfile.fetch(until_eof=True):
                rl = x.query_length
                ed = False
                for t in x.tags:
                    if t[0] == 'NM':
                        ed = t[1]
                        self.align_e(ed, rl, x)
                        break
    

        bamfile.close()
        self.outfile.close()
    

if __name__ == '__main__':

    descript = """
    
    Filter BAM/SAM alignment file by percent identity or edit distance.
    Output file format will always be sorted bam. 

    Incompatible flags:
    * -b and -s
    * -p and -g

    Examples:
    -b BAM.bam -o BAM.id97.bam (suggested option)
    -b BAM -p 0.97 -l 50 -t 5
    -s SAM -p 0.95 -l 150 -t 5
    -b BAM -g 1 3 -t 5 -m g
    
"""

    parser = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter, usage=argparse.SUPPRESS)
    parser.add_argument('-b', metavar='', type=str, nargs=1, default=[''], help='input BAM file (sorted or unsored)')
    parser.add_argument('-s', metavar='', type=str, nargs=1, default=[''], help='input SAM file')
    parser.add_argument('-o', metavar='', type=str, nargs=1, default=[''], help='output filtered file')
    parser.add_argument('-p', metavar='', type=float, nargs=1, default=[0.97], help='percent identity cutoff [0.97] (default)')
    parser.add_argument('-e', metavar='', type=str, nargs=1, default=[''], help='maximum edit distance (mismatch+gap+insert+delete)')
    parser.add_argument('-l', metavar='', type=int, nargs=1, default=[50], help='minimum length per read [50]')
    parser.add_argument('-i', metavar='', type=int, nargs=1, default=[''], help='file containing list of scaffolds of interest')
    parser.add_argument('-t', metavar='', type=int, nargs=1, default=[1], help='threads (SAM->BAM, BAM sorting) [1]')

    args = parser.parse_args()
    #
    bam = args.b[0]
    sam = args.s[0]
    out = args.o[0]
    interest = args.i[0]
    percent = args.p[0]
    edit = args.e[0]
    threads = args.t[0]
    if threads < 1:
        threads = 1
    length = args.l[0]

    if (sam and bam) or (not sam and not bam):
        print('\nError: provide SAM or BAM. Exiting.\n')
        exit()
    if sam and sam[-4:].lower() != '.sam':
        print('\nError: provide SAM must end in ".sam". Exiting.\n')
        exit()
    if bam and bam[-4:].lower() != '.bam':
        print('\nError: provide BAM must end in ".bam". Exiting.\n')
        exit()

    if edit:
        edit = int(edit)
        method = 'e'
    else:
        method = 'p'


    FilterCoverage(sam, bam, out, percent, edit, threads, length, method, interest)