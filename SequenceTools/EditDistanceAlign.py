from SequenceTools.Sequence import Sequence

class EditDistanceAlign:
    '''
    Uses edit distance to align two Sequence objects upon initialisation
        Attributes:
            sequenceA (Sequence)
            sequenceB (Sequence)
            dynamicProgramingMatrix (list)
    '''
    def __init__(self, sequenceA, sequenceB):
        '''
        Constructs the EditDistanceAlign instance using two sequences and aligns them

            Parameters:
                sequenceA (Sequence): The first sequence to align
                sequenceB (Sequence): The second sequence to align
        '''
        self.sequenceA = sequenceA
        self.sequenceB = sequenceB
        self.dynamicProgramingMatrix = []
        #We iterate over the two sequences as rows and columns in dynamicProgramingMatrix, constructing it as we go. The matrix is 1 field larger for pure indels at the beginning
        for i in range(0, self.sequenceA.getSequenceLength() + 1):
            self.dynamicProgramingMatrix.append([])
            for j in range(0, self.sequenceB.getSequenceLength() + 1):
                self.dynamicProgramingMatrix[i].append(0)
                #The first row and column are to be made purely of inserts / deletes, thus their values increase by 1
                if (i == 0):
                    self.dynamicProgramingMatrix[i][j] = j
                    continue
                if (j == 0):
                    self.dynamicProgramingMatrix[i][j] = i
                    continue
                #Upon match we take the diagonal
                if self.sequenceA.getSequence()[i - 1] == self.sequenceB.getSequence()[j - 1]:
                    self.dynamicProgramingMatrix[i][j] = self.dynamicProgramingMatrix[i - 1][j - 1]
                else:
                    #Upon mismatch we insert / delete / edit depending on which is less costly
                    self.dynamicProgramingMatrix[i][j] = min(self.dynamicProgramingMatrix[i - 1][j] + 1, self.dynamicProgramingMatrix[i - 1][j - 1] + 1, self.dynamicProgramingMatrix[i][j - 1] + 1)

    def getEditDistance(self):
        '''
        Returns the edit distance computed for this alignment

            Returns:
                distance (int): The edit distance of sequenceA and sequenceB
        '''
        #The distance is at the last column of the last row of dynamicProgramingMatrix
        return self.dynamicProgramingMatrix[self.sequenceA.getSequenceLength()][self.sequenceB.getSequenceLength()]

    def getAlignment(self):
        alignmentA = ''
        alignmentB = ''
        i = self.sequenceA.getSequenceLength()
        j = self.sequenceB.getSequenceLength()
        #We process the whole length of our sequences during backtracking to construct our alifned sequences
        while (i > 0) or (j > 0):
            #If we reach the end of one sequence, the rest are just indels while we reconstruct the rest of the other sequence
            if i == 0:
                alignmentA = '-' + alignmentA
                alignmentB = self.sequenceB.getSequence()[j - 1] + alignmentB
                j -= 1
                continue
            if j == 0:
                alignmentB = '-' + alignmentB
                alignmentA = self.sequenceA.getSequence()[i - 1] + alignmentA
                i -= 1
                continue
            #Upon match we backtrack diagonally
            if self.sequenceA.getSequence()[i - 1] == self.sequenceB.getSequence()[j - 1]:
                alignmentA = self.sequenceA.getSequence()[i - 1] + alignmentA
                alignmentB = self.sequenceB.getSequence()[j - 1] + alignmentB
                i -= 1
                j -= 1
                continue
            #Upon mismatch we backtrack to the lowest above value, continues used to avoid nesting elses
            if self.dynamicProgramingMatrix[i - 1][j - 1] < self.dynamicProgramingMatrix[i - 1][j] and self.dynamicProgramingMatrix[i - 1][j - 1] < self.dynamicProgramingMatrix[i][j - 1]:
                alignmentA = self.sequenceA.getSequence()[i - 1] + alignmentA
                alignmentB = self.sequenceB.getSequence()[j - 1] + alignmentB
                i -= 1
                j -= 1
                continue
            if self.dynamicProgramingMatrix[i - 1][j] < self.dynamicProgramingMatrix[i][j - 1]:
                alignmentB = '-' + alignmentB
                alignmentA = self.sequenceA.getSequence()[i - 1] + alignmentA
                i -= 1
                continue
            alignmentA = '-' + alignmentA
            alignmentB = self.sequenceB.getSequence()[j - 1] + alignmentB
            j -= 1
        return [alignmentA, alignmentB]
