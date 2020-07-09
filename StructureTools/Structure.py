from Bio import PDB
from StructureTools.Model import Model
from StructureTools.Atom import Atom
from StructureTools.Chain import Chain
from StructureTools.Residue import Residue
from StructureTools.Record import Record
from scipy.spatial import ConvexHull
from Exceptions.LoadException import LoadException
import math

class Structure(Record):
    '''
    A class representing a structure
        Attributes:
            record (Bio.PDB.Structure)
    '''
    def __init__(self, uri = None, id = None):
        '''
        Constructs the Structure object

            Parameters:
                uri (str): The string representation of the path to structure
                id (str): An optional id to represent the structure, otherwise it is taken from the filename sans extension
        '''
        self.record = None
        self.parent = None
        if not (uri is None):
            self.loadPdb(uri, id)

    def loadPdb(self, uri, id = None):
        '''
        Constructs the Structure from a PDB file

            Parameters:
                uri (str): The string representation of the path to structure
                id (str): An optional id to represent the structure, otherwise it is taken from the filename sans extension
        '''
        self.checkLoaded(False)
        if id is None:
            id = uri.split('.')[0].capitalize()
        parser = PDB.PDBParser(PERMISSIVE = True, QUIET = True)
        self.record = parser.get_structure(id, uri)

    def checkLoaded(self, shouldBeLoaded = True):
        '''Checks whether the object is loaded and raises exception when it is in the wrong state'''
        if self.isLoaded() != shouldBeLoaded:
            if shouldBeLoaded:
                raise LoadException('Structure is not loaded')
            else:
                raise LoadException('Structure already loaded')

    def isLoaded(self):
        return not (self.record is None)

    #Simple getters do not need comments
    def getStructure(self):
        self.checkLoaded()
        return self.record

    def getModelIds(self):
        self.checkLoaded()
        toReturn = []
        for model in self.record.get_models():
            toReturn.append(model.get_id())
        return toReturn

    def getModels(self):
        self.checkLoaded()
        toReturn = []
        for model in self.record.get_models():
            toReturn.append(Model(model, self))
        return toReturn

    def getChildren(self):
        self.checkLoaded()
        return self.getModels()

    def getModel(self, modelId):
        self.checkLoaded()
        for model in self.record.get_models():
            if model.get_id() == modelId:
                return Model(model, self)
        return None

    def getNumberOfModels(self):
        self.checkLoaded()
        return len(self.getModelIds())

    def getNumberOfChains(self):
        self.checkLoaded()
        count = 0
        for model in self.getModels():
            count += model.getNumberOfChains()
        return count

    def getNumberOfResidues(self):
        self.checkLoaded()
        count = 0
        for model in self.getModels():
            count += model.getNumberOfResidues()
        return count

    def getNumberOfAtoms(self):
        self.checkLoaded()
        count = 0
        for model in self.getModels():
            count += model.getNumberOfAtoms()
        return count

    def getAllAtoms(self):
        self.checkLoaded()
        atoms = []
        for model in self.getModels():
            for chain in model.getChains():
                for residue in chain.getResidues():
                    atoms += residue.getAtoms()
        return atoms

    def fastAtomsNearby(self, atom, distance):
        '''
        Uses the quick library function to find nearby atoms
            Parameters:
                atom (Atom): The atom around which to search
                distance (float): the distance in which to search
            Returns
                atoms (list)
        '''
        self.checkLoaded()
        searcher = PDB.NeighborSearch(list(self.record.get_atoms()))
        search = searcher.search(atom.getCoords(), distance, 'A')
        toReturn = []
        for atom in search:
            toReturn.append(Atom(atom, Residue(atom.get_parent(), Chain(atom.get_parent().get_parent(), Model(atom.get_parent().get_parent().get_parent(), self)))))
        return toReturn

    def fastResiduesNearby(self, atom, distance):
        '''
        Uses the quick library function to find nearby residues
            Parameters:
                atom (Atom): The atom around which to search
                distance (float): the distance in which to search
            Returns
                residues (list)
        '''
        self.checkLoaded()
        searcher = PDB.NeighborSearch(list(self.record.get_atoms()))
        search = searcher.search(atom.getCoords(), distance, 'R')
        toReturn = []
        for residue in search:
            toReturn.append(Residue(residue, Chain(residue.get_parent(), Model(residue.get_parent().get_parent(), self))))
        return toReturn

    #Could get even faster width in 2D using rotating calipers, but I have no idea as to what their equivalent is in 3D
    def fastWidth(self):
        '''
        Computes the width of the structure using scipy. First we get the convex hull of our atoms, upon which all our most distant candidates lie
            The fast proof is that if a candidate would be inside, there would be an even further one out on the hull
            Returns:
                width (float)
        '''
        self.checkLoaded()
        points = []
        for atom in self.getAllAtoms():
            points.append(atom.getCoords())
        hull = ConvexHull(points)

        def computeDistance(pointA, pointB):
            return math.sqrt(math.pow(pointA[0] - pointB[0], 2) + math.pow(pointA[1] - pointB[1], 2) + math.pow(pointA[2] - pointB[2], 2))

        maxDistance = 0
        for vertexA in hull.vertices:
            for vertexB in hull.vertices:
                distance = computeDistance(points[vertexA], points[vertexB])
                if  distance > maxDistance:
                    maxDistance = distance

        return maxDistance

    def getWidth(self):
        '''
        Computes the width of the structure hungrily by comparing all the atom distances
            Returns:
                width (float)
        '''
        self.checkLoaded()
        maxDistance = 0
        maxAtomA = None
        maxAtomB = None
        atoms = self.getAllAtoms()
        for atomA in atoms:
            for atomB in atoms:
                distance = atomA.getDistanceFrom(atomB)
                if distance > maxDistance:
                    maxDistance = distance
                    maxAtomA = atomA
                    maxAtomB = atomB
        return maxDistance
