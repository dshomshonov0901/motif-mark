# 1. FastaRecord class: Represents one FASTA sequence record, finds exon/intron segments based on case and reports length of sequence
# 2. Exon/IntronInternval: treats each exon and intron block as an interval with start and end coordinates
# 3. MotifPattern: holds the motif pattern that we check  against sequence window
# 4. MotifLocation: holds the locations of a motif found (good for storing coordinates for drawing)
# 5. MotifMark: actually runs the function and outputs image with Cairo

class FastaRecord:
    def __init__(self,record):
        self.record = record
    
    def __len__(self):
        for line in self.record.readlines():
            if not line.startswith(">"):
                return len(line)
    
    def __iter__(self):
        return iter(self.record)
    
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
    
    # def intron(self):
    #     down_intron = ""
    #     up_intron = ""
    #     for line in self.record.readlines():
    #         if not line.startswith(">"):
    #             for nuc in line:
    #                 if not nuc.isupper():
                        


record = open("test_fasta.txt", "r")

test = FastaRecord(record)
print(test.__len__())
print(test.exon_interval())


