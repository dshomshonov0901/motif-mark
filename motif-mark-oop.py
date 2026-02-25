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
            return self.header[1:].split()[0] #gene name
        return self.header.split()[0] if self.header else "record"

    def __len__(self):
        return len(self.sequence)

    def normalized_sequence(self):
        #uppercase all nucleotides so that motif finding is not case based
        return self.sequence.upper()

    def exon_intron_intervals(self):
        seq = self.sequence

        def feature_type(ch):
            #use in rendering
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
                #make sure the motif pattern is valid
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

#class MotifMark:
def __init__(self, fasta_path, motifs_path):
    self.fasta_path = fasta_path
    self.motifs_path = motifs_path

def run(self):
    records = self.read_fasta()
    motifs = self.read_motifs()
    hits_by_record = self.find_all_hits(records, motifs)
    self.render(records, motifs, hits_by_record)

def read_fasta(fasta_path):
    #reads in the fasta file and creates a list of FasaRecord objects per header in FASTA
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
    #reads in the motifs file and creates a list of MotifPattern objects per motif in the file
    motifs = []
    with open(motifs_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            motifs.append(MotifPattern(line))
    return motifs

def find_all_hits(records, motifs):
    #for FastaRecord object, looks for each motif hit and creates MotifLocation lists per hit
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
    #sorts hit starts numerically, assigns the MotifLocation lane parameter based on if the current hits start is greater than or equal
    #to the end of the lane
    #sort by start
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
            h.lane = len(lane_ends) #first lane is set
            lane_ends.append(h.end)
    return hits

def merge_hits_for_record(hits_for_one_record):
    """
    Merge overlapping/adjacent motif hits of the same motif on a single record so theres not too many lanes
    """
    record_name = hits_for_one_record[0].record_name
    grouped = {}
    for hit in hits_for_one_record:
        key = hit.motif.pattern
        if key not in grouped:
            grouped[key] = []
        grouped[key].append(hit)

    merged_hits = []

    for motif_key, hits in grouped.items():
        # Sort by start (then end)
        hits.sort(key=lambda h: (h.start, h.end))

        # Start with the first interval
        cur_start = hits[0].start
        cur_end = hits[0].end

        # Merge intervals that overlap or touch when next.start == cur_end, which we also merge to make one continuous bar
        for h in hits[1:]:
            if h.start <= cur_end:
                if h.end > cur_end:
                    cur_end = h.end
            else:
                # No overlap then finalize current interval, start a new one
                merged_hits.append(MotifLocation(hits[0].motif, record_name, cur_start, cur_end, 0))
                cur_start = h.start
                cur_end = h.end

        # Finalize the last merged interval for this motif
        merged_hits.append(MotifLocation(hits[0].motif, record_name, cur_start, cur_end, 0))

    # sort merged hits overall by start
    merged_hits.sort(key=lambda h: (h.start, h.end))

    return merged_hits


def render(records, hits_by_record, out_png, motifs):
        left_margin = 180
        right_margin = 30
        top_margin = 50
        drawable_width = 1000
        track_gap = 70
        height = top_margin + len(records) * track_gap + 120
        #making px_per_base to scale for the longest record in the fasta
        longest = max(len(r) for r in records)
        px_per_base = drawable_width / longest
        #painting white background
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,drawable_width + right_margin + left_margin, height)
        context = cairo.Context(surface)
        context.set_source_rgb(1,1,1)
        context.paint() #to make the whole background we start with paint not stroke

         # set colors for diff motifs, max 5 motifs so max 5 colors
        palette5 = [
            (0.18, 0.45, 0.80),
            (0.80, 0.25, 0.25),
            (0.20, 0.65, 0.35),
            (0.55, 0.35, 0.75),
            (0.90, 0.55, 0.15),
        ]
        motif_patterns = [m.pattern for m in motifs]
        if len(motif_patterns) > 5:
            raise ValueError("Max of 5 motifs supported.")
        
        motif_colors = {pat: palette5[i] for i, pat in enumerate(motif_patterns)}

        
        for row, record in enumerate(records):
            #adding gene names to left side
            baseline_y = top_margin + row * track_gap + 25
            x0 = left_margin
            x1 = left_margin + len(record) * px_per_base

            context.set_source_rgb(0,0,0)
            context.set_line_width(3)
            context.move_to(x0, baseline_y)
            context.line_to(x1, baseline_y)
            context.stroke()
            #choosing font
            context.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            context.set_font_size(12)
            rec_name = record.name()
            # label
            context.set_source_rgb(0, 0, 0)
            context.move_to(20, baseline_y + 4)
            context.show_text(rec_name)

            #exon rectangles
            intervals = record.exon_intron_intervals()

            for interval in intervals:
                if interval.feature_type == "exon":
                    x = left_margin + interval.start * px_per_base
                    w = (interval.end - interval.start) * px_per_base
                    exon_height = 18
                    y = baseline_y - exon_height / 2
                    h = exon_height
                    context.set_source_rgb(0, 0, 0)
                    context.rectangle(x, y, w, h)
                    context.fill()
                #drawing motifs
                hits = hits_by_record.get(record.name(), [])
                for hit in hits:
                    x = left_margin + hit.start * px_per_base
                    w = (hit.end - hit.start) * px_per_base
                    motif_height = 10
                    motif_gap = 3 #vertical gap between motifs in diff lanes
                    lane_offset = hit.lane * (motif_height + motif_gap)
                    motif_spacing_from_exon = 8
                    exon_height = 18
                    y = baseline_y - exon_height/2 - motif_spacing_from_exon - motif_height - lane_offset
                    r, g, b = motif_colors.get(hit.motif.pattern, (0.2, 0.2, 0.8))
                    context.set_source_rgb(r, g, b)
                    context.rectangle(x, y, w, motif_height)
                    context.fill()


        # legend
        legend_x = left_margin + 25
        legend_y = (height - 2 *top_margin)
        swatch = 12
        row_h = 18

        context.set_source_rgb(0, 0, 0)
        context.set_font_size(13)
        context.move_to(legend_x, legend_y)
        context.show_text("Motifs")
        context.set_font_size(12)

        y = legend_y + 18
        for pat in motif_patterns:
            r, g, b = motif_colors[pat]
            context.set_source_rgb(r, g, b)
            context.rectangle(legend_x, y - swatch + 2, swatch, swatch)
            context.fill()
            context.set_source_rgb(0, 0, 0)
            context.move_to(legend_x + swatch + 10, y + 2)
            context.show_text(pat)
            y += row_h

        surface.write_to_png(out_png)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", required=True, help="FASTA file")
    parser.add_argument("-m", required=True, help="motifs text file")
    args = parser.parse_args()

    fasta_path = args.f
    motifs_path = args.m

    #same fasta and png prefix 
    prefix, _ = os.path.splitext(os.path.basename(fasta_path))
    out_png = prefix + ".png"

    records = read_fasta(fasta_path)
    motifs = read_motifs(motifs_path)

    hits_by_record = find_all_hits(records, motifs)

    # merge + lanes per record
    for rec in records:
        name = rec.name()
        hits = hits_by_record.get(name, [])
        if hits:
            merged = merge_hits_for_record(hits)
            hits_by_record[name] = assign_lanes(merged)

    render(records, hits_by_record, out_png, motifs)

if __name__ == "__main__":
    main()
