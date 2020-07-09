class Record:
    '''
    A class representing an abstract record in the structure pipeline
        Attributes:
            record (object)
    '''
    def __init__(self, record, parent):
        self.record = record
        self.parent = parent

    def getId(self):
        return self.record.get_id()

    def __str__(self):
        return 'Record of type: ' + str(type(self)) + ' id: ' + str(self.getId())

    def getParent(self):
        return self.parent

    #Slow algos, but my own
    def getAtomsNearby(self, atom, distance):
        '''
        Traverses all the recorded atoms to find nearby atoms
            Parameters:
                atom (Atom): The atom around which to search
                distance (float): the distance in which to search
            Returns
                atoms (list)
        '''
        toReturn = []
        for child in self.getChildren():
            toReturn += child.getAtomsNearby(atom, distance)
        return toReturn

    def getResiduesNearby(self, atom, distance):
        '''
        Traverses all the recorded atoms to find nearby residues
            Parameters:
                atom (Atom): The atom around which to search
                distance (float): the distance in which to search
            Returns
                residues (list)
        '''
        toReturn = []
        for child in self.getChildren():
            toReturn += child.getResiduesNearby(atom, distance)
        return toReturn

    def getChildren(self):
        pass

    def getFullId(self):
        '''Returns path of selections to record'''
        return self.record.get_full_id()
