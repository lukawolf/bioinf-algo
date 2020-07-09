from Bio.Align import substitution_matrices

class ScoringMatrix:
    '''Represents a scoring matrix for MSA'''
    def __init__(self, uri = None, id = None):
        '''
        Constructs the ScoringMatrix object

            Parameters:
                uri (str): The string representation of the path to scoring matrix
                id (str): The matrix identifier
        '''
        self.record = None
        if not (uri is None):
            self.loadTable(uri, id)

    def loadTable(self, uri, id = None):
        '''
        Loads scoring matrix from file

            Parameters:
                uri (str): The string representation of the path to MSA
                id (str): The matrix identifier
        '''
        self.checkLoaded(False)
        if id is None:
            id = uri.split('.')[0].capitalize()
        self.id = id
        self.record = substitution_matrices.read(uri)

    #Simple getters do not need comments
    def getId(self):
        return self.id

    def getTable(self):
        return self.record

    def checkLoaded(self, shouldBeLoaded = True):
        '''Checks whether the object is loaded and raises exception when it is in the wrong state'''
        if self.isLoaded() != shouldBeLoaded:
            if shouldBeLoaded:
                raise LoadException('Structure is not loaded')
            else:
                raise LoadException('Structure already loaded')

    def isLoaded(self):
        return not (self.record is None)

    def score(self, a, b, universalPlaceholder = '*', universalPenalty = -10):
        '''
        Scores the two given characters according to the loaded scoring matrix
            Parameters:
                a (str): The first character to score
                b (str): The character to score against
                universalPlaceholder (str): The character the table's creators used to mark as universal match
                universalPenalty (int): The value to penalise when a or b is not covered by the table when there is no universalPlaceholder
        '''
        self.checkLoaded()
        a = a.upper()
        b = b.upper()
        line = self.record.get(a, None)
        if line is None:
            line = self.record.get(universalPlaceholder, None)
        if line is None:
            return universalPenalty
        toReturn = line.get(b, None)
        if toReturn is None:
            toReturn = line.get(universalPlaceholder, None)
        if toReturn is None:
            return universalPenalty
        return toReturn
