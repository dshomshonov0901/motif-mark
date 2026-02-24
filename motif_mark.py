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
    def __init__(self, header, sequence):
        self.header = header.strip()
        self.sequence = "".join(sequence.split())
        if not self.sequence:
            raise ValueError("Empty sequence for record: " + self.header)

    def name(self):
        if self.header.startswith(">"):
            return self.header[1:].split()[0]
        return self.header.split()[0] if self.header else "record"

    def __len__(self):
        return len(self.sequence)

    def normalized_sequence(self):
        return self.sequence.upper()

    def exon_intron_intervals(self):
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
    
class MotifPattern:
    """
    Holds a motif pattern and can test if a sequence window matches it,
    including IUPAC ambiguous nucleotide codes.
    """
    IUPAC = {
        "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
        "U": {"T"}, "R": {"A", "G"}, "Y": {"C", "T"}, "B": {"C", "G", "T"},
        "D": {"A", "G", "T"}, "K": {"G", "T"}, "M": {"A", "C"}, "H": {"A", "C", "T"},
        "V": {"A", "G", "C"}, "S": {"C", "G"}, "W": {"A", "T"}, "N": {"A", "G", "C", "T"}}
    
    def __init__(self, pattern):
            self.pattern = pattern.strip().upper()
            for ch in self.pattern:
                if ch not in self.IUPAC:
                    raise ValueError("Unsupported motif character: " + ch)

    def width(self):
        return len(self.pattern)
    
    def matches(self, window):
        window = window.upper()
        for ch1, ch2 in zip(self.pattern, window):
            if ch2 not in self.IUPAC[ch1]:
                return False
        return True
    
class MotifLocation:
    def __init__(self, motif, record_name, start, end, lane):
        self.motif = motif
        self.record_name = record_name
        self.start = start
        self.end = end
        self.lane = lane

    def length(self) -> int:
        return self.end - self.start

class MotifMark:
    def __init__(self, fasta_path, motifs_path):
        self.fasta_path = fasta_path
        self.motifs_path = motifs_path

    def run(self):
        records = self.read_fasta()
        motifs = self.read_motifs()
        hits_by_record = self.find_all_hits(records, motifs)
        self.render(records, motifs, hits_by_record)

    def read_fasta(fasta_path):
        records = []
        current_header = None
        current_sequence = []
        with open(fasta_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header:
                        full_sequence = "".join(current_sequence)
                        records.append(FastaRecord(current_header, full_sequence))
                    current_header = line
                    current_sequence = []
                else:
                    current_sequence.append(line)
        records.append(FastaRecord(current_header, "".join(current_sequence)))
        return records

    def read_motifs(motifs_path):
        motifs = []
        with open(motifs_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                motifs.append(MotifPattern(line))
        return motifs

    def find_all_hits(records, motifs):
        hits_by_record = {}
        for r in records:
            record_hits = []
            seq = r.normalized_sequence()
            for m in motifs:
                width = m.width()
                for i in range(0,len(seq) - width +1):
                    window = seq[i:i+width]
                    if m.matches(window) == True:
                        record_hits.append(MotifLocation(m,r,i,i+width,0))
            hits_by_record[r.name()] = record_hits
        return hits_by_record

    def assign_lanes(hits_for_one_record):
        hits = sorted(hits_for_one_record, key=lambda h: h.start)
        lane_ends = []
        for h in hits:
            placed = False

            for i in range(len(lane_ends)):
                if h.start >= lane_ends[i]:
                    h.lane = i
                    lane_ends[i] = h.end
                    placed = True
                    break

            if not placed:
                h.lane = len(lane_ends)
                lane_ends.append(h.end)

        return hits

