#Todo: think about disordered atom class

from StructureTools.Record import Record
import math

class Atom(Record):
    '''
    A class representing an atom
        Attributes:
            record (Bio.PDB.Atom)
    '''
    #Simple getters do not need comments
    def getChildren(self):
        return None

    def getCoords(self):
        return self.record.get_coord()

    def getDistanceFrom(self, atom):
        '''Computes basic 3D distance between self and another atom instance'''
        myCoords = self.getCoords()
        theirCoords = atom.getCoords()
        return math.sqrt(math.pow(myCoords[0] - theirCoords[0], 2) + math.pow(myCoords[1] - theirCoords[1], 2) + math.pow(myCoords[2] - theirCoords[2], 2))

    def getAtomsNearby(self, atom, distance):
        #As the final link in getting nearby atoms, it can only return self and that only when it is not too far
        toReturn = []
        if self.getDistanceFrom(atom) <= distance:
            toReturn.append(self)
        return toReturn

    def isDisordered(self):
        return self.record.is_disordered()

    def getSerialNumber(self):
        return self.record.get_serial_number()

    def getName(self):
        return self.record.get_name()

    def getBFactor(self):
        return self.record.get_bfactor()

    def getOccupancy(self):
        return self.record.get_occupancy()
