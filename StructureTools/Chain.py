from StructureTools.Residue import Residue
from StructureTools.Record import Record
from StructureTools.Atom import Atom
from SequenceTools.Sequence import Sequence
from StructureTools.Residue import Residue
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Chain(Record):
    '''
    A class representing a chain
        Attributes:
            record (Bio.PDB.Chain)
    '''
    #Simple getters do not need comments
    def getChain(self):
        return self.record

    def getResidues(self):
        toReturn = []
        for residue in self.record.get_residues():
            toReturn.append(Residue(residue, self))
        return toReturn

    def getChildren(self):
        return self.getResidues()

    def getResidueIds(self):
        toReturn = []
        for residue in self.record.get_residues():
            toReturn.append(residue.get_id()[1])
        return toReturn

    def getResidue(self, residueId):
        for residue in self.record.get_residues():
            if residue.get_id()[1] == residueId:
                return Residue(residue, self)
        return None

    def getNumberOfResidues(self):
        return len(self.getResidueIds())

    def getAtomsOfType(self, type):
        toReturn = []
        for atom in self.record.get_atoms():
            if atom.get_id() == type:
                toReturn.append(Atom(atom, Residue(atom.get_parent(), self)))
        return toReturn

    def getNumberOfAtoms(self):
        count = 0
        for residue in self.getResidues():
            count += residue.getNumberOfAtoms()
        return count

    def toSequence(self):
        residueString = ""
        for residue in self.getResidues():
            residueString += residue.getSingleLetterName()
        record = SeqRecord(Seq(residueString))
        sequence = Sequence()
        sequence.loadRecord(record)
        #If the first ID is above 1 then the creator is aware of the issue, else we compute our own padding
        sequence.needsPadding = (self.getResidueIds()[0] <= 1)
        return sequence
