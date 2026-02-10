1. FastaRecord class: Represents one FASTA sequence record, finds exon/intron segments based on case and reports length of sequence
2. Exon/IntronInternval: treats each exon and intron block as an interval with start and end coordinates
3. MotifPattern: holds the motif pattern that we check  against sequence window
4. MotifLocation: holds the locations of a motif found (good for storing coordinates for drawing)
5. MotifMark: actually runs the function and outputs image with Cairo

