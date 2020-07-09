from Bio import SeqIO
from Exceptions.LoadException import LoadException

class Sequence:
    '''
    Represents a sequence
        Attributes:
            record (SeqRecord)
            needsPadding (bool): Normally sequence needs padding when used from structure and is compared to MSA. Flag set by constructor in chain if needed
    '''
    def __init__(self, uri = None):
        '''
        Constructs the Sequence object

            Parameters:
                uri (str): The string representation of the path to sequence
        '''
        self.record = None
        self.needsPadding = False
        if not (uri is None):
            self.loadFasta(uri)

    def loadFasta(self, uri, recordIndex = None):
        '''
        Loads sequence from file

            Parameters:
                uri (str): The string representation of the path to sequence
                recordIndex (int): The sequence index within the fasta file to consider if there are more records (if None the user will be asked in console to specify)
        '''
        self.checkLoaded(False)
        records = list(SeqIO.parse(uri, "fasta"))
        selectedIndex = 0
        self.record = records[selectedIndex]
        if (len(records) > 1):
            if (recordIndex is None):
                while True:
                    try:
                        selectedIndex = int(input("There are " + str(len(records)) + " records in this fasta file. Which record do you want to load?\n"))
                        self.record = records[selectedIndex - 1]
                    except:
                        print('Please enter a valid number between 1 and ' + str(len(records)) + "\n")
                        continue
                    break
            else:
                self.record = records[recordIndex]

    def loadRecord(self, record):
        '''Loads the sequence from a given BioPython record'''
        self.checkLoaded(False)
        self.record = record

    def checkLoaded(self, shouldBeLoaded = True):
        '''Checks whether the object is loaded and raises exception when it is in the wrong state'''
        if self.isLoaded() != shouldBeLoaded:
            if shouldBeLoaded:
                raise LoadException('Sequence is not loaded')
            else:
                raise LoadException('Sequence already loaded')

    #Simple getters do not need comments
    def isLoaded(self):
        return not (self.record is None)

    def getId(self):
        self.checkLoaded()
        return self.record.id

    def getName(self):
        self.checkLoaded()
        return self.record.name

    def getDescription(self):
        self.checkLoaded()
        return self.record.description

    def getFeatures(self):
        self.checkLoaded()
        return self.record.features

    def getSequence(self):
        self.checkLoaded()
        return self.record.seq

    def getSequenceLength(self):
        self.checkLoaded()
        return len(self.record.seq)

    def getSubSequence(self, fromChar, toChar):
        '''Gets a string representation of a subsequence within this sequence'''
        self.checkLoaded()
        if fromChar < 0 or toChar >= self.getSequenceLength() or fromChar > toChar:
            raise IndexError('From or to is out of range')
        return self.getSequence()[fromChar:toChar + 1]

    def matchChainToMsaSequences(self, msa):
        '''Matches a chain created sequence to a msa sequence requiring exact match except for gaps and tolerating shorter structure sequences via padding'''
        matches = []
        for msaSequence in msa.getAllSequences():
            noSpaceMsa    = str(msaSequence.getSequence()).replace("-", "")
            paddingNeeded = len(noSpaceMsa) - len(self.getSequence())
            for paddingBefore in range(0, paddingNeeded + 1):
                paddingAfter = paddingBefore - paddingNeeded
                compareTo = None
                if paddingAfter < 0:
                    compareTo = noSpaceMsa[paddingBefore:paddingAfter]
                else:
                    compareTo = noSpaceMsa[paddingBefore:]
                if self.getSequence() == compareTo:
                    spacesBeforeId = []
                    spacesFound = 0
                    for i in range(0, msaSequence.getSequenceLength()):
                        if msaSequence.getSequence()[i] == "-":
                            spacesFound += 1
                        else:
                            spacesBeforeId.append(spacesFound)
                    matches.append({'paddingBefore': paddingBefore, 'paddingAfter': paddingAfter, 'spacesBeforeId': spacesBeforeId, 'matchedMsaSequence': msaSequence})
        return matches


    def mapChainIdToMsaPos(self, id, match):
        #If id is that large, it is out of bounds
        if id >= len(match['spacesBeforeId']):
            return None
        #-1 as sequences are 1 indexed and we use 0 indexed operations
        return id + match['spacesBeforeId'][id] + (int(self.needsPadding) * match['paddingBefore']) - 1
