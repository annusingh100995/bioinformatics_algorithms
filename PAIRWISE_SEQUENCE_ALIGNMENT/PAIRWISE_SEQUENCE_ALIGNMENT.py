import numpy as np 
from MySeq import MySeq

# Objective function 
# Function to create a substitution matrix based on the values of
# a match or mismatch score


class SequenceAlignment():
    
    def __init__(self,sm,g):
        # gap panelties
        self.g = g
        #substitution matrix
        self.sm = sm
        # the score matrix from alignment
        self.S = None
        # matrix to backtrack the sequence
        self.T = None
        self.seq1 = None
        self.seq2 = None
        
    def score_position(self,c1,c2):
        if c1 == "_" or c2 == "_":
            return self.g
        else:
            return self.sm[c1+c2]
    
    def score_align(self,seq1,seq2):
        res = 0
        for i in range(len(seq1)):
            res += self.score_position(seq1[i], seq2[i])
        return res
    
    def needleman_wunsch(self,seq1,seq2):
        self.S = [[0]]
        self.T = [[0]]
        
        for col in range(1,len(seq2)+1):
            self.S[0].append(self.g*col)
            self.T[0].append(3)
        
        for row in range(1,len(seq1)+1):
            self.S.append([self.g*row])
            self.T.append([2])
        
        for row in range(1,len(seq1)+1):
            for col in range(1,len(seq2)+1):
                s1 = self.S[row-1][col-1] + self.score_position(seq1[row-1],seq2[col-1])
                s2 = self.S[row][col-1] + self.g
                s3 = self.S[row-1][col] + self.g
                self.S[row].append(max(s1,s2,s3))
                #self.T[row].append(self.max_index(s1,s2,s3))
                self.T[row].append(np.argmax([s1,s2,s3]))
        return (self.S, self.T)
    
    def max_index(self,v1,v2,v3):
        if v1 > v2:
            if v1 > v3:
                return 1
            else:
                return 3
        else:
            if v2 > v3 :
                return 2
            else:
                return 3
    
    def recover_alignment(self,seq1,seq2):
        result = ["",""]
        i = len(seq1)-1
        j = len(seq2)-1
        while i>0 or j>0:
            if self.T[i][j] == 1:
                
                result[0] = seq1[i-1]+result[0]
                result[1] = seq2[j-1]+result[1]
                i -= 1
                j -= 1
            elif self.T[i][j] == 3:
                
                result[0] = "_" + result[0]
                result[1] = seq2[j-1] + result[1]
                j -= 1
            else:
                
                result[0] = seq1[i-1] + result[0]
                result[1] = "_" + result[1]
                i -= 1
            
        return result
    
    def local_alignment(self,seq1,seq2):
        self.S = [[0]]
        self.T = [[0]]
        max_score = 0
        for col in range(1,len(seq1)+1):
            self.S[0].append(col*self.g)
            self.T[0].append(0)
        
        for row in range(1,len(seq1)+1):
            self.S.append([0])
            self.T.append([0])
        
        for row in range(1,len(seq1)+1):
            for col in range(1,len(seq2)+1):
                s1 = self.S[row-1][col-1]+ self.score_position(seq1[row-1],seq2[col-1])
                s2 = self.S[row][col-1]+self.g
                s3 = self.S[row-1][col-1]+self.g
                b = max(s1,s2,s3)
                if b <= 0:
                    self.S[row].append(0)
                    self.T[row].append(0)
                else:
                    self.S[row].append(b)
                    self.T[row].append(np.argmax([s1,s2,s3]))
                    if b > max_score:
                        max_score = b
        
        return (self.S, self.T, max_score)
    
    def recover_align_local(self,seq1,seq2):
        result = ["LOCAL1","LOCAL2"]
        row , col = self.max_mat(self.S)
        while self.T[row][col]>0:
            if self.T[row][col] == 1:
                result[0] = seq1[row-1] + result[0]
                result[1] = seq2[col-1] + result[1]
                row -= 1
                col -= 1
            elif self.T[row][col] == 3:
                result[0] = "_" + result[0]
                result[1] = seq2[col-1] + result[1]
                j -= 1
            elif self.T[row][col] == 2:
                result[0] = seq1[row-1] + result[0]
                result[1] = "_" + result[1]
                row -= 1
        return result
    
    def max_mat(self,mat):
        maxval = mat[0][0]
        maxrow = 0
        maxcol = 0
        for i in range(0,len(mat)):
            for j in range(0,len(mat[i])):
                if mat[i][j] > maxval:
                    maxval = mat[i][j]
                    maxrow = i
                    maxcol = j
        return (maxrow,maxcol)
    
def create_subs_matrix(match,mismatch,alphabet):
    sm = {}
    for c1 in alphabet:
        for c2 in alphabet:
            if(c1==c2):
                sm[c1+c2] = match
            else:
                sm[c1+c2] = mismatch
    return sm

def tedt_DNA():
    sm= create_subs_matrix(1,0,'ATGC')
    print(sm)


tedt_DNA()


def test():
    seq1 = "AAAATATGCAAA"
    seq2 = "AAAGAAATGCAA"
    #seq1 = "_CAGTGCATG_ACATA"
    #seq2 = "TCAG_GC_TCTACAGA"
    sm = create_subs_matrix(6,-2,"ACGT")
    g = -3
    seq_align = SequenceAlignment(sm,-8)
    mat = seq_align.score_position(seq1[0],seq2[0])
    print(mat)
    res = seq_align.needleman_wunsch(seq1,seq2)
    S = res[0]
    T = res[1]
    #print(S)
    print("T")
    #print(T)
    align  = seq_align.recover_alignment(seq1,seq2)
    print(align[0])
    print(align[1])
    print("LOCAL RESULTS")
    local_align = seq_align.recover_align_local(seq1,seq2)
    print(S)
    print(T)
    print(local_align[0])
    print(local_align[1])

test()