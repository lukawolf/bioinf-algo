from StructureTools.Chain import Chain
from StructureTools.Record import Record

class Model(Record):
    '''
    A class representing a model
        Attributes:
            record (Bio.PDB.Model)
    '''
    #Simple getters do not need comments
    def getModel():
        return self.record

    def getChains(self):
        toReturn = []
        for chain in self.record.get_chains():
            toReturn.append(Chain(chain, self))
        return toReturn

    def getChildren(self):
        return self.getChains()

    def getChainIds(self):
        toReturn = []
        for chain in self.record.get_chains():
            toReturn.append(chain.get_id())
        return toReturn

    def getChain(self, chainId):
        for chain in self.record.get_chains():
            if chain.get_id() == chainId:
                return Chain(chain, self)
        return None

    def getNumberOfChains(self):
        return len(self.getChainIds())

    def getNumberOfResidues(self):
        count = 0
        for chain in self.getChains():
            count += chain.getNumberOfResidues()
        return count

    def getNumberOfAtoms(self):
        count = 0
        for chain in self.getChains():
            count += chain.getNumberOfAtoms()
        return count
