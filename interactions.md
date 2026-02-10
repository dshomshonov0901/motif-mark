1. **MotifMark** reads the FASTA file and creates one **FastaRecord** object per sequence entry
2. **MotifMark** reads the motifs file and creates one **MotifPattern** object per motif line.
3. For each **FastaRecord** and each **MotifPattern**, **MotifMark** scans along the sequence. At every position it extracts a window and calls **MotifPattern** to test for a match
4. When a match is found, **MotifMark** creates a **MotifLocation** object storing the motif identity and the start/end coordinates on that **FastaRecord**
5. **MotifMark** uses pycairo to draw one figure containing all sequence tracks. For each **FastaRecord**, it draws introns/exons using the **ExonIntronInterval** coordinates, then draws motif marks using the **MotifLocation** coordinates
