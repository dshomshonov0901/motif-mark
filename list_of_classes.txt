FastaRecord class: Represents one FASTA sequence record, finds exon/intron segments based on case and reports length of sequence
Exon/IntronInternval: treats each exon and intron block as an interval with start and end coordinates
MotifPattern: holds the motif pattern that we check  against sequence window
MotifLocation: holds the locations of a motif found (good for storing coordinates for drawing)
MotifMark: actually runs the function and outputs image with Cairo

