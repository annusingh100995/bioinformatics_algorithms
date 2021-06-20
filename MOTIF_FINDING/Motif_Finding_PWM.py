from random import randint
from random import random

def create_matrix_book(nrows,ncols):
    res = [ ]
    for i in range (0, nrows):
        res.append([0]*(ncols))
    return res

def print_matrix_book(mat):
    for i in range (0, len (mat)): print (mat[i])

# Class for biological sequences

class MySeq:
    
    def __init__(self,seq,seq_type='DNA'):
        self.seq = seq.upper()
        self.seq_type = seq_type
    
    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self,n):
        return self.seq[n]
    
    def __getslice__(self,i,j):
        return self.seq[i:j]
    
    def __str__(self):
        return self.seq
    
    def get_seq_biotype(self):
        return self.seq_type
    
    def show_info_seq(self):
        print("Sequence:",self.seq, "biotype:",self.seq_type)
    
    def alphabet(self):
        if(self.seq_type == 'DNA'):
            return 'ACGT'
        elif(self.seq_type == "RNA"):
            return 'ACGU'
        elif(self.seq_type == 'PROTEIN'):
            return 'ACDEFGHIKLMNPQRSTVWY'
        else:
            return None
    
    def validate(self):
        alp = self.alphabet()
        res = True
        i = 0
        while i <len(self.seq) and res:
            if self.seq[i] not in alp :
                res = False
            else:
                i += 1
        return res
    
    def transcription(self):
        if (self.seq_type == 'DNA'):
            return MySeq(self.seq.replace('T','U'),"RNA")
        else:
            return None
    
    def reverse_comp(self):
        if (self.seq_type != 'DNA'):
            return None
        comp = ""
        for c in self.seq:
            if (c == 'A'): comp = 'T'+comp
            elif (c == 'T'): comp = 'A'+comp
            elif (c == 'G'): comp = 'C'+comp
            elif (c == 'C'): comp = 'G'+comp
        return MySeq(comp, "DNA")
    
    def translate(self, iniPos = 0):
        if (self.seq_type != 'DNA'):
            return None
        seq_aa = ''
        for pos in range(iniPos,len(self.seq)-2, 3):
            cod = self.seq[pos:pos+3]
            seq_aa += translate_codon(cod)
        return MySeq(seq_aa, "PROTEIN")
    
    def translate_codon(self,cod):
       
        """Translates a codon into an aminoacid using an internal
        dictionary with the standard genetic code."""
        
        tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", "TGT":"C", "TGC":"C","GAT":"D", "GAC":"D",
        "GAA":"E", "GAG":"E","TTT":"F", "TTC":"F","GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G","CAT":"H", "CAC":"H",
        "ATA":"I", "ATT":"I", "ATC":"I","AAA":"K", "AAG":"K",
        "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L","ATG":"M", "AAT":"N", "AAC":"N",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T","GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V","TGG":"W",
        "TAT":"Y", "TAC":"Y","TAA":"_", "TAG":"_", "TGA":"_"}
        
        if cod in tc:
            return tc[cod]
        else :
            return None
    
    def translate_seq(dna_seq, ini_pos = 0):
        """ Translates a DNA sequence into an aminoacid sequence. """
        #"Invalid DNA sequence"
        
        assert dna_seq.validate()
        
        seqm = dna_seq
        
        seq_aa = ""
        
        for pos in range(ini_pos, len(seqm)-2, 3):
            cod = seqm[pos:pos+3]
            seq_aa += seqm.translate_codon(cod)
        
        return seq_aa


def test_MySeq():
    s = MySeq("ATGTGATAAGAATAGAATGCTGAATAAATAGAATGACAT")
    s1 = MySeq("MKVVLSVQERSVVSLL", "PROTEIN")
    print("Test Slicing",s[2:5])
    print('The Sequence :', s.show_info_seq())
    print("Protein Seq:",s1.show_info_seq())
    print("Validating ",s.validate())
    print('Reverse Compliment:',s.reverse_comp())
    print('Transcription:',s.transcription())
    print("translation:", s.translate_seq())
    print("DNA SEQ Alphabets:",s.alphabet())
    print("PROTEIN SEQ Alphabets:",s1.alphabet())

test_MySeq()

# Creating motif class
class MyMotifs:
    
    def __init__(self,seqs = [],pwm = [],symbol = None):
        
        
        if seqs:
            # Length of the sequences,sequenced are of equal length
            #self.size = len(seqs[0])
            self.size = len(seqs[0])
            # the input sequences
            self.seqs = seqs
            # The alphabets of the sequence
            self.symbol = seqs[0].alphabet()
            # frequency count mattrix for the symbols for seqs
            self.do_counts()
            # position weight matrices of the seqs
            self.create_pwm()
        else:
            self.pwm = pwm
            self.size = len(pwm[0])
            self.symbol = symbol
        
    def __len__(self):
        return self.size
    
    # Function to count the frequency for each symbol for the motif
    def do_counts(self):
        
        self.counts = create_matrix_book(len(self.symbol), self.size)
        for s in self.seqs:
            for i in range(0,self.size):
                lin = self.symbol.index(s[i])
                self.counts[lin][i] = self.counts[lin][i]+1
                
    # Func. to create position weight matrix
    def create_pwm(self):
        if self.counts == None:
            self.do_counts()
        self.pwm = create_matrix_book(len(self.symbol), self.size)
        for i in range(len(self.symbol)):
            for j in range(self.size):
                self.pwm[i][j] = float(self.counts[i][j]/len(self.seqs[0]))
    
    """This function scans every position of the motif and
    returns the most frequent symbol at each position."""
    def consensus(self):
        res = ''
        for col in range(self.size):
            maxcol = self.counts[0][col]
            maxcol_index = 0
            for row in range(1,len(self.symbol)):
                if self.counts[row][col] > maxcol:
                    maxcol= self.counts[row][col]
                    maxcol_index = row
            res += self.symbol[maxcol_index]
        return res
    
    """This is similar to concensus only difference being,it returns 
    the most frequent symbol only if it is at least 50% of the input sequence."""
    
    def masked_consensus(self):
        res = ''
        for col in range(self.size):
            maxcol = self.counts[0][col]
            maxcol_index = 0
            
            for row in range(1,len(self.symbol)):
                if self.counts[row][col] > maxcol:
                    maxcol = self.counts[row][col]
                    maxcol_index = row
                if maxcol > len(self.seqs)/2:
                    res += self.symbol[maxcol_index]
                else:
                    res += '_'
        return res
    
    # Probability for a particular sequence
    def probability_sequence(self,seq):
        res = 1.0
        for col in range(self.size):
            lin = self.symbol.index(seq[col])
            res *= self.pwm[lin][col]
        return res
    
    # Find probabilty of all the sub sequences possibel in a seq
    def probability_all_positions(self,seq):
        res = []
        for k in range(len(seq)-self.size):
            res.append(self.probability_sequence(seq[k:k+self.size]))
        return res
    
    # Return the index of the most probable sub-sequence of the seq
    def most_probable_sequence(self,seq):
        maximum = -1.0
        max_index = -1
        for k in range(len(seq) - self.size):
            p = self.probability_sequence(seq[k:k+self.size])
            if(p > maximum):
                maximum = p
                max_index = k
        return max_index
    
    def create_motif(self,seq):
        #from Myseq import Myseq
        l = []
        for s in seqs:
            ind = self.most_probable_sequence(s.seq)
            subseq = Myseq(s[ind:self.size],s.get_seq_biotype())
            l.append(subseq)
        return MyMotifs(l)


def test():
    seq1 = MySeq("AAAGTT")
    seq2 = MySeq("CACGTG")
    seq3 = MySeq("TTGGGT")
    seq4 = MySeq("GAAAGT")
    seq5 = MySeq("AACCAT")
    seq6 = MySeq("AACCCT")
    seq7 = MySeq("AAACCT")
    seq8 = MySeq("GAACCT")
    lseqs = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8]
    motifs = MyMotifs(lseqs)
    print ("Counts matrix")
    print_matrix_book(motifs.counts)
    print ("PWM")
    print_matrix_book(motifs.pwm)
    print ("Sequence alphabet")
    print (motifs.symbol)
    [ print (s) for s in lseqs]
    print ("Consensus sequence")
    print (motifs.consensus())
    print ("Masked Consensus sequence")
    print (motifs.masked_consensus())
    print ("Prob of seq : AAACCT: ",motifs.probability_sequence("AAACCT"))
    print ("Prob of seq : ATACAG: ",motifs.probability_sequence("ATACAG"))
    print (motifs.most_probable_sequence("CTATAAACCTTACATC"))
    print("Motif size = ", motifs.size)
    print("prob all pos output:",motifs.probability_all_positions(seq="CTATAAACCTTACATC"))


test()

class StochasticsMotifFinding:
    
    def __init__(self,size = 4,seqs = None):
        self.motif_size = size
        if(seqs != None):
            self.seqs = seqs
            self.alphabet = seqs[0].alphabet()
        else:
            self.seqs = seqs
    
    def __len__(self):
        return len(self.seqs)
    
    def __getitem__(self,n):
        return self.seqs[n]
    
    def seq_size(self,i):
        return len(self.seqs[i])
    
    def create_motif_from_indicies(self,index_vector):
        pseqs =[]
        for i,ind in enumerate(index_vector):
            pseqs.append(MySeq(self.seqs[i][ind:(ind+self.motif_size)],self.seqs[i].get_seq_biotype()))
        return MyMotifs(pseqs)
    
    
    # Scoring function
    def score_add(self,s):
        score = 0
        motif = self.create_motif_from_indicies(s)
        mat = motif.counts
        for col in range(len(mat[0])):
            maxcol = mat[0][col]
            for row in range(len(mat)):
                if mat[row][col] > maxcol:
                    maxcol = mat[row][col]
                score += maxcol
        return score
    
    def score_multiplicative(self,s):
        score = 1.0
        motif = self.create_motif_from_indicies(s)
        mat = motif.counts
        for col in range(len(mat[0])):
            maxcol = mat[0][col]
            for row in range(len(mat)):
                if mat[row][col] > maxcol:
                    maxcol = mat[row][col]
                score *= maxcol
        return score
    
    def create_random_index(self):
        from random import randint
        
        s = [0]*len(self.seqs)
        for k in range(len(self.seqs)):
            s[k] = randint(0,self.seq_size(k)-self.motif_size)
        
        return s
        
    def heuristics_stochastics(self):
        #from random import randint
        
        #s = [0]*len(self.seqs)
        #for k in range(len(self.seqs)):
        #    s[k] = randint(0,self.seq_size(k)-self.motif_size)
        
        s = self.create_random_index()
        
        motif = self.create_motif_from_indicies(s)
        motif.create_pwm()
        sc = self.score_multiplicative(s)
        bestsol = s
        improve = True
        while(improve):
            for k in range(len(s)):
                s[k] = motif.most_probable_sequence(self.seqs[k])
            
            if self.score_add(s) > sc:
                sc = self.score_add(s)
                bestsol = s
                motif = self.create_motif_from_indicies(s)
                motif.create_pwm()
            else:
                improve = False
            return bestsol
        
    def gibbsamp_motif_finding(self, iteration = 100):
        s = self.create_random_index()
        best_initial_position = list(s)
        best_score = self.score_multiplicative(s)
        for i in range(iteration):
            # random index to remove the corresponding seq
            seq_index = randint(0,len(self.seqs)-1)
            seq_removed = self.seqs[seq_index]
            # removing the initial pos.of the removed index
            s.pop(seq_index)
            # these are all the seq after the selected seq is removed
            removed_seq = self.seqs.pop(seq_index)
            motif = self.create_motif_from_indicies(s)
            motif.create_pwm()
            # Inserting the removed seq back to the original pos.
            self.seqs.insert(seq_index,seq_removed)
            # prob. of all the motif possible in the removeed seq
            r = motif.probability_all_positions(seq_removed)
            # select a random but also highly prob motif from the removed seq
            pos = self.roulette(r)
            # inserting the new pos removed from the initial pos. back
            s.insert(seq_index, pos)
            score = self.score_multiplicative(s)
            if score > best_score:
                best_score = score
                best_initial_position = list(s)
        return best_initial_position
    
    """ This function is used to create a sample space
    for the experiment of spinning once a roulette wheel.
    All the positions have a chance of being picked,
    but this chance is proportional to their score"""
    def roulette(self,f):
        from random import random
        total = 0.0
        for x in f:
            total += (0.01+x)
        val = random()*total
        acum = 0.0
        idx = 0
        while acum < val:
            acum += (f[idx]+0.01)
            idx += 1
        return idx-1

def test_STOCHASTIC_MOTIF_FINDINF():
    seq1 = MySeq("ATAGAGCTGA","DNA")
    seq2 = MySeq("ACGTAGATGA","DNA")
    seq3 = MySeq("AAGATAGGGG","DNA")
    dna_seq = [seq1,seq2,seq3]
    em = StochasticsMotifFinding(size=5, seqs=dna_seq)
    print("Random Index")
    a = em.create_random_index()
    print(a)
    print("Additive score for the above index:",em.score_add(a))
    print("Multiplicative score for the above index:",em.score_multiplicative(a))
    print("Finding the best solution")
    s = em.heuristics_stochastics()
    best_motif = em.create_motif_from_indicies(s)
    print("The consencusmotif : ",best_motif.consensus())
    print("The masked consensus motif : ",best_motif.masked_consensus())
    print("Best intial pos. for the motif in each seq: ", em.gibbsamp_motif_finding())

def test2_STOCHASTIC_MOTIF_FINDINF():
    seq1 = MySeq("CTGCGGAGTTCGCGAC","DNA")
    seq2 = MySeq("GAATGAACACACTGCG","DNA")
    seq3 = MySeq("ATGCGATTCAGCATGT","DNA")
    seq4 = MySeq("TTCGAAGGTAGCGATG","DNA")
    seq5 = MySeq("CGTAGGCGTCGGCACT","DNA")
    seq6 = MySeq("CCGGGATGTAAAGTCA","DNA")
    seq7 = MySeq("ACTTGTTTCTCTTATG","DNA")
    seq8 = MySeq("ATCCACCGTCGGCACT","DNA")
    seq9 = MySeq("GTACTTTTCGGCTACG","DNA")
    seq_list = [seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8,seq9]
    em = StochasticsMotifFinding(size=5, seqs=seq_list)
    print("Random Index")
    a = em.create_random_index()
    print(a)
    print("Additive score for the above index:",em.score_add(a))
    print("Multiplicative score for the above index:",em.score_multiplicative(a))
    print("Finding the best solution")
    s = em.heuristics_stochastics()
    best_motif = em.create_motif_from_indicies(s)
    print("The consencus motif : ",best_motif.consensus())
    print("The masked consensus motif : ",best_motif.masked_consensus())
    print("Best intial pos. for the motif in each seq: ", em.gibbsamp_motif_finding())

for i in range(10):
    print(test2_STOCHASTIC_MOTIF_FINDINF())
