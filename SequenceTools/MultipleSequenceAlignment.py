from Bio import AlignIO
from Exceptions.LoadException import LoadException
from SequenceTools.Sequence import Sequence
from sortedcontainers import SortedList

class MultipleSequenceAlignment:
    '''
    Represents multiple sequence alignment (MSA)
        Attributes:
            record (MultipleSeqAlignment)
            id (string)
    '''
    def __init__(self, uri = None):
        '''
        Constructs the MultipleSequenceAlignment object

            Parameters:
                uri (str): The string representation of the path to MSA
        '''
        self.record = None
        self.id = None
        if not (uri is None):
            self.loadClustal(uri)

    def isLoaded(self):
        '''
        Checks whether the MSA is loaded

            Returns:
                isLoaded(bool)
        '''
        return not (self.record is None)

    def checkLoaded(self, shouldBeLoaded = True):
        '''Checks whether the object is loaded and raises exception when it is in the wrong state'''
        if self.isLoaded() != shouldBeLoaded:
            if shouldBeLoaded:
                raise LoadException('MSA is not loaded')
            else:
                raise LoadException('MSA already loaded')

    def loadClustal(self, uri):
        '''
        Loads clustal format MSA from file

            Parameters:
                uri (str): The string representation of the path to MSA
        '''
        self.checkLoaded(False)
        self.id = uri.split('.')[0].capitalize()
        self.record = AlignIO.read(uri, "clustal")

    #Simple getters do not need comments
    def getId(self):
        return self.id

    def getMsa(self):
        self.checkLoaded()
        return self.record

    def getNumberOfSequences(self):
        return len(self.record)

    def getSequenceByPosition(self, position):
        '''
            Gets sequence according to its position in MultipleSequenceAlignment

                Parameters:
                    position (int)

                Returns:
                    sequence (Sequence)
        '''
        self.checkLoaded()
        if (position < 0) or (position > (self.getNumberOfSequences() - 1)):
            raise IndexError('Position ' + str(position) + ' is out of range: 0 to ' + str(self.getNumberOfSequences() - 1))
        sequence = Sequence()
        sequence.loadRecord(self.record[position])
        return sequence

    def getSequenceIds(self):
        '''Gets the list of all sequence available IDs'''
        self.checkLoaded()
        toReturn = []
        for sequence in self.record:
            toReturn.append(sequence.id)
        return toReturn

    def getSequenceById(self, id):
        '''
            Gets sequence according to its ID in MultipleSequenceAlignment

                Parameters:
                    id (str)

                Returns:
                    sequence (Sequence)
        '''
        self.checkLoaded()
        for sequence in self.record:
            if sequence.id == id:
                toReturn = Sequence()
                toReturn.loadRecord(sequence)
                return toReturn
        return None

    def getNumberOfColumns(self):
        self.checkLoaded()
        return self.record.get_alignment_length()

    def getAllSequences(self):
        '''
            Gets all sequences in MultipleSequenceAlignment

                Returns:
                    sequences (list)
        '''
        self.checkLoaded()
        toReturn = []
        for sequence in self.record:
            toReturnSequence = Sequence()
            toReturnSequence.loadRecord(sequence)
            toReturn.append(toReturnSequence)
        return toReturn

    def getColumn(self, position):
        '''
        Gets the residues in a given MSA column
            Parameters:
                position (int): Position of the column
            Returns:
                sequences (list)
        '''
        self.checkLoaded()
        if (position < 0) or (position >= self.getNumberOfColumns()):
            raise IndexError('Position ' + str(position) + ' is out of range: 0 to ' + str(self.getNumberOfColumns() - 1))
        toReturn = []
        for sequence in self.record:
            toReturn.append(sequence[position])
        return toReturn

    def getColumnSumOfPairs(self, position, scoringMatrix):
        '''
        Gets the SP score in a given MSA column
            Parameters:
                position (int): Position of the column
                scoringMatrix (ScoringMatrix)
            Returns:
                score (int)
        '''
        self.checkLoaded()
        column = self.getColumn(position)
        score = 0
        for x in range(0, len(column)):
            for y in range(x + 1, len(column)):
                score += scoringMatrix.score(column[x], column[y])
        return score

    def getColumnsSumOfPairs(self, scoringMatrix):
        '''
        Gets the SP score in all MSA columns
            Parameters:
                scoringMatrix (ScoringMatrix)
            Returns:
                scores (list)
        '''
        self.checkLoaded()
        score = []
        for i in range(0, self.getNumberOfColumns()):
            score.append(self.getColumnSumOfPairs(i, scoringMatrix))
        return score

    def getMsaSumOfPairs(self, scoringMatrix):
        '''
        Gets the sum of the SP scores in all MSA columns
            Parameters:
                scoringMatrix (ScoringMatrix)
            Returns:
                score (int)
        '''
        self.checkLoaded()
        score = 0
        for i in range(0, self.getNumberOfColumns()):
            score += self.getColumnSumOfPairs(i, scoringMatrix)
        return score

    def getConservationScoreForSequence(self, sequenceId, scoringMatrix):
        '''
        Gets the SP score for a MSA sequence ignoring the spaces in it
            Parameters:
                sequenceId (object): int or str, uses getSequenceByPosition or getSequenceById to acquire the sequence
                scoringMatrix (ScoringMatrix)
            Returns:
                score (int)
        '''
        self.checkLoaded()
        sequence = None
        if isinstance(sequenceId, int):
            sequence = self.getSequenceByPosition(sequenceId)
        else:
            sequence = self.getSequenceById(sequenceId)
        score = 0
        for i in range(0, sequence.getSequenceLength()):
            if sequence.getSubSequence(i, i) != '-':
                score += self.getColumnSumOfPairs(i, scoringMatrix)
        return score

    def getConservationScoresForSubSequence(self, fromPosition, toPosition, scoringMatrix):
        '''
        Gets the SP scores for MSA columns in between and including given positions
            Parameters:
                fromPosition (int)
                toPosition (int)
                scoringMatrix (ScoringMatrix)
            Returns:
                scores (list)
        '''
        self.checkLoaded()
        if (fromPosition < 0) or (fromPosition > (self.getNumberOfColumns() - 1)) or (toPosition < 0) or (toPosition > (self.getNumberOfColumns() - 1)) or (fromPosition > toPosition):
            raise IndexError('Position is out of range')
        toReturn = []
        for i in range(fromPosition, toPosition + 1):
            toReturn.append(self.getColumnSumOfPairs(i, scoringMatrix))
        return toReturn

    def getTopConservedSpans(self, amount, scoringMatrix):
        '''
        Gets the top spans in the MSA by sum of SP scores
            Parameters:
                amount (int)
                scoringMatrix (ScoringMatrix)
            Returns:
                scores (list): List of dictionaries containing score, positionFrom and positionTo
        '''
        self.checkLoaded()
        #Scores are kept in ordered lists for ease of keeping the correct amount
        scores = SortedList(key = lambda entry: entry['score'])
        columnScores = self.getColumnsSumOfPairs(scoringMatrix)
        #Checks all ranges and sums the score for them. A sliding window would also be nice, but that would require a fixed length of spans or iterating over the lengths
        for positionFrom in range(0, self.getNumberOfColumns()):
            for positionTo in range(positionFrom, self.getNumberOfColumns()):
                score = 0
                for i in range(positionFrom, positionTo + 1):
                    score += columnScores[i]
                scores.add({'score': score, 'positionFrom': positionFrom, 'positionTo': positionTo})
                if len(scores) > amount:
                    scores.remove(scores[0])
        return scores

    def getTopConservedResidues(self, amount, scoringMatrix):
        '''
        Gets the top positions in the MSA by SP scores
            Parameters:
                amount (int)
                scoringMatrix (ScoringMatrix)
            Returns:
                scores (list): List of dictionaries containing score and position
        '''
        self.checkLoaded()
        #Scores are kept in ordered lists for ease of keeping the correct amount
        scores = SortedList(key = lambda entry: entry['score'])
        columnScores = self.getColumnsSumOfPairs(scoringMatrix)
        for i in range(0, self.getNumberOfColumns()):
            scores.add({'score': columnScores[i], 'position': i})
            if len(scores) > amount:
                scores.remove(scores[0])
        return scores
