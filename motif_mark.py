#!/usr/bin/env python3

# 1. FastaRecord class: Represents one FASTA sequence record, finds exon/intron segments based on case and reports length of sequence
# 2. Exon/IntronInternval: treats each exon and intron block as an interval with start and end coordinates
# 3. MotifPattern: holds the motif pattern that we check  against sequence window
# 4. MotifLocation: holds the locations of a motif found (good for storing coordinates for drawing)
# 5. MotifMark: actually runs the function and outputs image with Cairo

import argparse
import os
import cairo

class ExonIntronInterval:
    """
    Represents one contiguous exon or intron block with start and end corrdinates and feature type
    """
    def __init__(self, start, end, feature_type):
        self.start = start
        self.end = end
        self.feature_type = feature_type

    def length(self):
        return self.end - self.start

    def __repr__(self):
        return f"{self.feature_type}({self.start}, {self.end})"

class FastaRecord:
    """
    Represents one FASTA record, header + sequence. Exons are uppercase bases. Introns are lowercase bases.
    """
    def __init__(self,header, sequence):
        self.header = header.strip()
        self.sequence = "".join(sequence.split()) #removes newlines and whitespace
    
    def name(self):
        # label for plot, gets gene name
        if self.header:
            return self.header.split()[0]
        return "record"
    
    def __len__(self):
        return len(self.sequence)
    
    def __iter__(self):
        return iter(self.record)
    
    def normalized_sequence(self):
        # for motif scanning we ignore exon/intron case and upper everything
        return self.sequence.upper()
    
    def exon_interval(self):
        exon = []
        self.record.seek(0)
        for line in self.record.readlines():
            for pos, nuc in enumerate(line):
                if not line.startswith(">"):
                    print(pos, line)
                    if nuc.isupper():
                        exon.append(pos)
        return exon
    
    def exon_intron_intervals(self):
        """
        Returns a list of ExonIntronInterval objects covering the whole sequence.
        Splits whenever case changes.
        """
        seq = self.sequence

        def feature_type(ch):
            return "exon" if ch.isupper() else "intron"

        intervals = []
        current_type = feature_type(seq[0])
        start = 0

        for i in range(1, len(seq)):
            t = feature_type(seq[i])
            if t != current_type:
                intervals.append(ExonIntronInterval(start, i, current_type))
                start = i
                current_type = t

        intervals.append(ExonIntronInterval(start, len(seq), current_type))
        return intervals
    
#testingggggg
if __name__ == "__main__":
    r = FastaRecord(">INSR chr19:7150261-7150808 (reverse complement)", "atgtccacatgtagtcacgtttgacatcccagggccacctcagcaggccgtctctggggagaattttctctgatttcttccccttcccttgctggacccctgcacctgctggggaagatgtagctcactccgtctagcaagtgatgggagcgagtggtccagggtcaaagccagggtgcccttactcggacacatgtggcctccaagtgtcagagcccagtggtctgtctaatgaagttccctctgtcctcaaaggcgttggttttgtttccacagAAAAACCTCTTCAGGCACTGGTGCCGAGGACCCTAGgtatgactcacctgtgcgacccctggtgcctgctccgcgcagggccggcggcgtgccaggcagatgcctcggagaacccaggggtttctgtggctttttgcatgcggcgggcagctgtgctggagagcagatgcttcaccaattcagaaatccaatgccttcactctgaaatgaaatctgggcatgaatgtggggagaaaccttcactaacacactcttgctaaaacatagaatca")
    print("name:", r.name())
    print("length:", len(r))
    print("sequence:", r.sequence)
    print("intervals:", r.exon_intron_intervals())

# record = open("test.fasta")

# test = FastaRecord(record)
# print(test.exon_interval())


