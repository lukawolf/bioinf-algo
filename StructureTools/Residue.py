#Todo: think about disordered residue class

from StructureTools.Atom import Atom
from StructureTools.Record import Record

mapping = {
    "ala":"A",
    "arg":"R",
    "asn":"N",
    "asp":"D",
    "asx":"B",
    "cys":"C",
    "glu":"E",
    "gln":"Q",
    "glx":"Z",
    "gly":"G",
    "his":"H",
    "ile":"I",
    "leu":"L",
    "lys":"K",
    "met":"M",
    "phe":"F",
    "pro":"P",
    "ser":"S",
    "thr":"T",
    "trp":"W",
    "tyr":"Y",
    "val":"V",
}

class Residue(Record):
    '''
    A class representing a residue
        Attributes:
            record (Bio.PDB.Residue)
    '''
    #Simple getters do not need comments
    def getId(self):
        return self.record.get_id()[1]

    def getSegmentId(self):
        return self.record.get_segid()

    def getAtoms(self):
        toReturn = []
        for atom in self.record.get_atoms():
            toReturn.append(Atom(atom, self))
        return toReturn

    def getChildren(self):
        return self.getAtoms()

    def getAtomIds(self):
        toReturn = []
        for atom in self.record.get_atoms():
            toReturn.append(atom.get_id())
        return toReturn

    def getAtom(self, atomId):
        for atom in self.record.get_atoms():
            if atom.get_id() == atomId:
                return Atom(atom, self)
        return None

    def getNumberOfAtoms(self):
        return len(self.getAtomIds())

    def getResiduesNearby(self, atom, distance):
        #As the final link in getting nearby residues it checks whether an atom in it is close enough and then returns the properly wrapped self
        for myAtom in self.getAtoms():
            if atom.getDistanceFrom(myAtom) <= distance:
                return [self]
        return []

    def isDisordered(self):
        return self.record.is_disordered()

    def getResidueName(self):
        return self.record.get_resname()

    def getSingleLetterName(self):
        return mapping.get(self.getResidueName().lower()) or ''
