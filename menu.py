from PyInquirer import style_from_dict, Token, prompt, Validator, ValidationError
from SequenceTools.Sequence import Sequence
from SequenceTools.hammingDistance import hammingDistance
from SequenceTools.EditDistanceAlign import EditDistanceAlign
from StructureTools.Structure import Structure
from SequenceTools.ScoringMatrix import ScoringMatrix
from SequenceTools.MultipleSequenceAlignment import MultipleSequenceAlignment
import os
from pathlib import Path

def cls():
    '''Clears the console'''
    if os.name == 'nt':
        _ = os.system('cls')
    else:
        _ = os.system('clear')

#Storage for menu loaded objects
sequences = {}
structures = {}
multipleSequenceAlignments = {}
scoringMatrices = {}

#A basic style for our menu
style = style_from_dict({
    Token.QuestionMark: '#E91E63 bold',
    Token.Selected: '#673AB7 bold',
    Token.Instruction: '',  # default
    Token.Answer: '#2196f3 bold',
    Token.Question: '',
})

class NumberValidator(Validator):
    '''Validates integers'''
    def validate(self, document):
        try:
            int(document.text)
        except ValueError:
            raise ValidationError(
                message='Please enter a valid number',
                cursor_position=len(document.text)
            )

class FloatValidator(Validator):
    '''Validates floats'''
    def validate(self, document):
        try:
            float(document.text)
        except ValueError:
            raise ValidationError(
                message='Please enter a valid number',
                cursor_position=len(document.text)
            )

class FileValidator(Validator):
    '''Validates file uris'''
    def validate(self, document):
        my_file = Path(document.text)
        if not my_file.is_file():
            raise ValidationError(
                message='Please enter a valid file uri',
                cursor_position=len(document.text)
            )

def doMenuAction(action):
    '''Holds and executes menu actions based on selections, its contents are mostly self-explainatory, but some subfunctions do get their own doc strings'''
    def invalid():
        print('Invalid action!')

    def exit():
        print('Leaving program')

    def loadTestingFiles():
        for uri in ['TestingFiles/3ijf.pdb', 'TestingFiles/3f3c.pdb', 'TestingFiles/3owe.pdb', ]:
            structure = Structure(uri)
            structures[structure.getId()] = structure
        for uri in ['TestingFiles/A9WZ33.fasta', 'TestingFiles/A9WZ33DIF.fasta', 'TestingFiles/A9WZ33Double.fasta',
        'TestingFiles/click.fasta', 'TestingFiles/clock.fasta', 'TestingFiles/lacks.fasta', ]:
            sequence = Sequence(uri)
            sequences[sequence.getId()] = sequence
        for uri in ['TestingFiles/Blosum62.txt', ]:
            sm = ScoringMatrix(uri)
            scoringMatrices[sm.getId()] = sm
        for uri in ['TestingFiles/CytidineDeaminase.clustal', 'TestingFiles/p53.clustal', ]:
            msa = MultipleSequenceAlignment(uri)
            multipleSequenceAlignments[msa.getId()] = msa

    def selectFile():
        subMenu = [
            {
                'type': 'input',
                'name': 'uri',
                'message': 'Enter file uri',
                'validate': FileValidator,
            },
        ]
        return prompt(subMenu, style=style)['uri']

    def loadFasta():
        sequence = Sequence(selectFile())
        sequences[sequence.getId()] = sequence

    def selectSequence():
        if len(sequences) < 1 :
            print('No sequences are loaded yet!')
            input('Press Enter to continue')
            return None
        subMenu = [
            {
                'type': 'list',
                'name': 'choice',
                'message': 'Which sequence?',
                'choices': sequences.keys(),
            }
        ]
        return prompt(subMenu, style=style)['choice']

    def getDescription():
        key = selectSequence()
        if key is None:
            return
        print(sequences[key].getDescription())
        input('Press Enter to continue')

    def getSequence():
        key = selectSequence()
        if key is None:
            return
        print(sequences[key].getSequence())
        input('Press Enter to continue')

    def getSequenceLength():
        key = selectSequence()
        if key is None:
            return
        print(sequences[key].getSequenceLength())
        input('Press Enter to continue')

    def getSubSequence():
        key = selectSequence()
        if key is None:
            return
        subMenu = [
            {
                'type': 'input',
                'name': 'fromChar',
                'message': 'From position 1 to ' + str(sequences[key].getSequenceLength()),
                'validate': NumberValidator,
                'filter': lambda val: int(val)
            },
            {
                'type': 'input',
                'name': 'toChar',
                'message': 'To position',
                'validate': NumberValidator,
                'filter': lambda val: int(val)
            },
        ]
        answers = prompt(subMenu, style=style)
        print(sequences[key].getSubSequence(answers['fromChar'] - 1, answers['toChar'] - 1))
        input('Press Enter to continue')

    def getHammingDistance():
        key1 = selectSequence()
        if key1 is None:
            return
        key2 = selectSequence()
        try:
            print(hammingDistance(sequences[key1], sequences[key2]))
        except ValueError as e:
            print(e)
        input('Press Enter to continue')

    def getEditDistanceAlign():
        key1 = selectSequence()
        if key1 is None:
            return
        key2 = selectSequence()
        alignment = EditDistanceAlign(sequences[key1], sequences[key2])
        print("The distance is: " + str(alignment.getEditDistance()))
        alignedSequences = alignment.getAlignment()
        print("The alignment is: \n" + alignedSequences[0] + "\n" + alignedSequences[1] + "\n")
        input('Press Enter to continue')

    def loadPdb():
        structure = Structure(selectFile())
        structures[structure.getId()] = structure

    def selectStructure():
        if len(structures) < 1 :
            print('No structures are loaded yet!')
            input('Press Enter to continue')
            return None
        subMenu = [
            {
                'type': 'list',
                'name': 'choice',
                'message': 'Which structure?',
                'choices': structures.keys(),
            }
        ]
        return prompt(subMenu, style=style)['choice']

    def structureDetails():
        key = selectStructure()
        if key is None:
            return None
        print("Structure ID: " + str(structures[key].getId()) + "\n")
        print("#models: " + str(structures[key].getNumberOfModels()))
        print("#chains: " + str(structures[key].getNumberOfChains()))
        print("#residues: " + str(structures[key].getNumberOfResidues()))
        print("#atoms: " + str(structures[key].getNumberOfAtoms()) + "\n")
        print("Structure width: " + str(structures[key].fastWidth()) + "\n")
        input('Press Enter to continue')

    def selectFromSet(set, name):
        if (set is None) or (len(set) < 1):
            print("No " + name + "s exist in selected path")
            input('Press Enter to continue')
            return None
        choices = []
        for element in set:
            choices.append({
                'name': str(element),
                'value': element,
            })
        subMenu = [
            {
                'type': 'list',
                'name': 'choice',
                'message': 'Which ' + name + '?',
                'choices': choices,
            }
        ]
        return prompt(subMenu, style=style)['choice']

    def selectModel():
        structureKey = selectStructure()
        if structureKey is None:
            return None
        modelKey = selectFromSet(structures[structureKey].getModelIds(), 'model')
        if modelKey is None:
            return None
        return structures[structureKey].getModel(modelKey)

    def modelDetails():
        model = selectModel()
        if model is None:
            return
        print("Model ID: " + str(model.getId()) + "\n")
        print("#chains: " + str(model.getNumberOfChains()))
        print("#residues: " + str(model.getNumberOfResidues()))
        print("#atoms: " + str(model.getNumberOfAtoms()) + "\n")
        input('Press Enter to continue')

    def selectChain():
        model = selectModel()
        if model is None:
            return None
        chainKey = selectFromSet(model.getChainIds(), 'chain')
        if chainKey is None:
            return None
        return model.getChain(chainKey)

    def chainDetails():
        chain = selectChain()
        if chain is None:
            return None
        print("Chain ID: " + str(chain.getId()) + "\n")
        print("#residues: " + str(chain.getNumberOfResidues()))
        print("#atoms: " + str(chain.getNumberOfAtoms()) + "\n")
        input('Press Enter to continue')

    def selectResidue():
        chain = selectChain()
        if chain is None:
            return None
        residueKey = selectFromSet(chain.getResidueIds(), 'residue')
        if residueKey is None:
            return None
        return chain.getResidue(residueKey)

    def residueDetails(residue = None, enterToContinue = True):
        if residue is None:
            residue = selectResidue()
        if residue is None:
            return None
        print("Residue ID: " + str(residue.getId()) + "\n")
        print("#atoms: " + str(residue.getNumberOfAtoms()) + "\n")
        print("Is disordered: " + str(residue.isDisordered()) + "\n")
        if enterToContinue:
            input('Press Enter to continue')

    def selectAtom():
        residue = selectResidue()
        if residue is None:
            return None
        atomKey = selectFromSet(residue.getAtomIds(), 'atom')
        if atomKey is None:
            return None
        return residue.getAtom(atomKey)

    def atomDetails(atom = None, enterToContinue = True):
        if atom is None:
            atom = selectAtom()
        if atom is None:
            return None
        print("Atom ID: " + str(atom.getId()) + "\n")
        print("Name: " + str(atom.getName()))
        print("Coordinates: X:" + str(atom.getCoords()[0]) + " Y: " + str(atom.getCoords()[1]) + " Z: " + str(atom.getCoords()[2]))
        print("B factor: " + str(atom.getBFactor()))
        print("Occupancy: " + str(atom.getOccupancy()))
        print("Is disordered: " + str(atom.isDisordered()) + "\n")
        if enterToContinue:
            input('Press Enter to continue')

    def atomsWithinDistance():
        atom = selectAtom()
        if atom is None:
            return None
        subMenu = [
            {
                'type': 'input',
                'name': 'distance',
                'message': 'Within distance: ',
                'validate': FloatValidator,
                'filter': lambda val: float(val)
            },
        ]
        distance = prompt(subMenu, style=style)['distance']
        structure = atom.getParent().getParent().getParent().getParent()
        atoms = structure.fastAtomsNearby(atom, distance)
        for atom in atoms:
            atomDetails(atom, False)
        input('Press Enter to continue')

    def residuesWithinDistance():
        atom = selectAtom()
        if atom is None:
            return None
        subMenu = [
            {
                'type': 'input',
                'name': 'distance',
                'message': 'Within distance: ',
                'validate': FloatValidator,
                'filter': lambda val: float(val)
            },
        ]
        distance = prompt(subMenu, style=style)['distance']
        structure = atom.getParent().getParent().getParent().getParent()
        residues = structure.fastResiduesNearby(atom, distance)
        for residue in residues:
            residueDetails(residue, False)
        input('Press Enter to continue')

    def loadMsa():
        msa = MultipleSequenceAlignment(selectFile())
        multipleSequenceAlignments[msa.getId()] = msa

    def selectMsa():
        if len(multipleSequenceAlignments) < 1 :
            print('No MSAs are loaded yet!')
            input('Press Enter to continue')
            return None
        subMenu = [
            {
                'type': 'list',
                'name': 'choice',
                'message': 'Which MSA?',
                'choices': multipleSequenceAlignments.keys(),
            }
        ]
        return prompt(subMenu, style=style)['choice']

    def getMsaSequenceCount():
        key = selectMsa()
        if key is None:
            return None
        print(multipleSequenceAlignments[key].getNumberOfSequences())
        input('Press Enter to continue')

    def getMsaLength():
        key = selectMsa()
        if key is None:
            return None
        print(multipleSequenceAlignments[key].getNumberOfColumns())
        input('Press Enter to continue')

    def getMsaSequenceIds():
        key = selectMsa()
        if key is None:
            return None
        for sequenceId in multipleSequenceAlignments[key].getSequenceIds():
            print(sequenceId)
        input('Press Enter to continue')

    def getMsaSequence():
        key = selectMsa()
        if key is None:
            return None
        sequenceKey = selectFromSet(multipleSequenceAlignments[key].getSequenceIds(), 'sequence')
        if sequenceKey is None:
            return None
        sequences[sequenceKey] = multipleSequenceAlignments[key].getSequenceById(sequenceKey)
        print('Sequence ' + sequenceKey + ' loaded into sequences')
        input('Press Enter to continue')

    def selectMsaColumn(key):
        subMenu = [
            {
                'type': 'input',
                'name': 'column',
                'message': 'Column number (1 through ' + str(multipleSequenceAlignments[key].getNumberOfColumns()) + '): ',
                'validate': NumberValidator,
                'filter': lambda val: int(val)
            },
        ]
        return prompt(subMenu, style=style)['column']

    def getMsaColumn():
        key = selectMsa()
        if key is None:
            return None
        column = selectMsaColumn(key)
        for character in multipleSequenceAlignments[key].getColumn(column - 1):
            print(character)
        input('Press Enter to continue')

    def loadScoringMatrix():
        sm = ScoringMatrix(selectFile())
        scoringMatrices[sm.getId()] = sm

    def selectScoringMatrix():
        if len(scoringMatrices) < 1 :
            print('No scoring matrices are loaded yet!')
            input('Press Enter to continue')
            return None
        subMenu = [
            {
                'type': 'list',
                'name': 'choice',
                'message': 'Which scoring matrix?',
                'choices': scoringMatrices.keys(),
            }
        ]
        return prompt(subMenu, style=style)['choice']

    def getColumnSumOfPairs():
        key = selectMsa()
        if key is None:
            return None
        smKey = selectScoringMatrix()
        if smKey is None:
            return None
        column = selectMsaColumn(key)
        print(multipleSequenceAlignments[key].getColumnSumOfPairs(column, scoringMatrices[smKey]))
        input('Press Enter to continue')

    def getMsaSumOfPairs():
        key = selectMsa()
        if key is None:
            return None
        smKey = selectScoringMatrix()
        if smKey is None:
            return None
        print(multipleSequenceAlignments[key].getMsaSumOfPairs(scoringMatrices[smKey]))
        input('Press Enter to continue')

    def getConservationScoreForSequence():
        key = selectMsa()
        if key is None:
            return None
        smKey = selectScoringMatrix()
        if smKey is None:
            return None
        sequenceKey = selectFromSet(multipleSequenceAlignments[key].getSequenceIds(), 'sequence')
        if sequenceKey is None:
            return None
        print(multipleSequenceAlignments[key].getConservationScoreForSequence(sequenceKey, scoringMatrices[smKey]))
        input('Press Enter to continue')

    def getConservationScoresForSubSequence():
        key = selectMsa()
        if key is None:
            return None
        smKey = selectScoringMatrix()
        if smKey is None:
            return None
        columnA = selectMsaColumn(key)
        columnB = selectMsaColumn(key)
        for score in multipleSequenceAlignments[key].getConservationScoresForSubSequence(columnA, columnB, scoringMatrices[smKey]):
            print(score)
        input('Press Enter to continue')

    def getTopConservedSpans():
        key = selectMsa()
        if key is None:
            return None
        smKey = selectScoringMatrix()
        if smKey is None:
            return None
        subMenu = [
            {
                'type': 'input',
                'name': 'count',
                'message': 'How many?',
                'validate': NumberValidator,
                'filter': lambda val: int(val)
            },
        ]
        count = prompt(subMenu, style=style)['count']
        for conservedSpan in multipleSequenceAlignments[key].getTopConservedSpans(count, scoringMatrices[smKey]):
            print("Conserved span from position: " + str(conservedSpan['positionFrom']) + " to position: " + str(conservedSpan['positionTo']) + " with score: " + str(conservedSpan['score']))
        input('Press Enter to continue')

    def getTopConservedResidues():
        key = selectMsa()
        if key is None:
            return None
        smKey = selectScoringMatrix()
        if smKey is None:
            return None
        subMenu = [
            {
                'type': 'input',
                'name': 'count',
                'message': 'How many?',
                'validate': NumberValidator,
                'filter': lambda val: int(val)
            },
        ]
        count = prompt(subMenu, style=style)['count']
        for conservedResidue in multipleSequenceAlignments[key].getTopConservedResidues(count, scoringMatrices[smKey]):
            print("Conserved residue at position: " + str(conservedResidue['position']) + " with score: " + str(conservedResidue['score']))
        input('Press Enter to continue')

    def compareConservation():
        '''Uses MultipleSequenceAlignment, ScoringMatrix and Structure to compare the active site conservation compared to other positions'''
        #Selecting each needed datapoint is self-explainatory
        structureKey = selectStructure()
        if structureKey is None:
            return None
        msaKey = selectMsa()
        if msaKey is None:
            return None
        smKey = selectScoringMatrix()
        if smKey is None:
            return None
        structure = structures[structureKey]
        model=structure.getModel(0)
        msa = multipleSequenceAlignments[msaKey]
        sm = scoringMatrices[smKey]
        chainKey = selectFromSet(model.getChainIds(), "chain")
        if chainKey is None:
            return None
        chain = model.getChain(chainKey)
        chainSequence = chain.toSequence()
        matches = chainSequence.matchChainToMsaSequences(msa)
        if len(matches) < 1:
            print('Structure does not match any MSA sequence')
            input('Press Enter to continue')
            return None
        subMenu = [
            {
                'type': 'input',
                'name': 'atomType',
                'message': 'Enter atom type (e.q. ZN)',
            },
        ]
        atomType = prompt(subMenu, style=style)['atomType']
        foundAtoms = chain.getAtomsOfType(atomType)
        if len(foundAtoms) < 1:
            print ('No such atoms were found')
            input('Press Enter to continue')
            return None
        foundAtom = None
        if len(foundAtoms) > 1:
            for found in foundAtoms:
                print("Atom coordinates: " + found.getCoords())
            subMenu = [
                {
                    'type': 'input',
                    'name': 'atom',
                    'message': 'Multiple atoms found, select atom index from 1 to ' + str(len(foundAtoms)),
                    'validate': NumberValidator,
                },
            ]
            foundAtom = foundAtoms[int(prompt(subMenu, style=style)['atom']) - 1]
        else:
            foundAtom = foundAtoms[0]
        subMenu = [
            {
                'type': 'input',
                'name': 'distance',
                'message': 'How distant residues are to be considered?',
                'validate': FloatValidator,
            },
        ]
        distance = float(prompt(subMenu, style=style)['distance'])
        #Here we finally can get nearby residues to our selected active site atom
        residuesNearby = structure.fastResiduesNearby(foundAtom, distance)
        residues = []
        for residue in residuesNearby:
            #We ignore water and out of bounds residues (HETATM)
            if (residue.getResidueName() != "HOH") and (residue.getParent().getId() == chainKey) and not (chainSequence.mapChainIdToMsaPos(residue.getId(), matches[0]) is None):
                residues.append(residue)
        if len(residues) < 1:
            print("No residues found within distance from selected atom")
            input('Press Enter to continue')
            return None
        scores = msa.getColumnsSumOfPairs(sm)
        #For our residues (marked by zero based ids) we check their scores against the scores of other MSA columns
        for residue in residues:
            residueMsaPosition = chainSequence.mapChainIdToMsaPos(residue.getId(), matches[0])
            if residueMsaPosition is None:
                continue
            print("Residue ID " + str(residue.getId()) + " MSA position: " + str(residueMsaPosition + 1) + " residue name: " + residue.getResidueName())
            print("MSA column: " + str(msa.getColumn(residueMsaPosition)))
            greater = 0
            for i in range(0, len(scores)):
                if scores[i] > scores[residueMsaPosition]:
                    print("\tGreater score position " + str(i + 1) + " (Difference from conserved score: " + str(scores[i] - scores[residueMsaPosition]) + ")")
                    print("\t\tMSA column: " + str(msa.getColumn(i)))
                    greater += 1
            print("There are " + str(greater) + " higher scores than active site residue\n")
        input('Press Enter to continue')

    #The dictionary of menu actions to their servicing functions
    switcher = {
        'loadTestingFiles': loadTestingFiles,
        'loadFasta': loadFasta,
        'exit': exit,
        'getDescription': getDescription,
        'getSequence': getSequence,
        'getSequenceLength': getSequenceLength,
        'getSubSequence': getSubSequence,
        'getHammingDistance': getHammingDistance,
        'getEditDistanceAlign': getEditDistanceAlign,
        'loadPdb': loadPdb,
        'structureDetails': structureDetails,
        'modelDetails': modelDetails,
        'chainDetails': chainDetails,
        'residueDetails': residueDetails,
        'atomDetails': atomDetails,
        'atomsWithinDistance': atomsWithinDistance,
        'residuesWithinDistance': residuesWithinDistance,
        'loadMsa': loadMsa,
        'getMsaSequenceCount': getMsaSequenceCount,
        'getMsaLength': getMsaLength,
        'getMsaSequenceIds': getMsaSequenceIds,
        'getMsaSequence': getMsaSequence,
        'getMsaColumn': getMsaColumn,
        'loadScoringMatrix': loadScoringMatrix,
        'getColumnSumOfPairs': getColumnSumOfPairs,
        'getMsaSumOfPairs': getMsaSumOfPairs,
        'getConservationScoreForSequence': getConservationScoreForSequence,
        'getConservationScoresForSubSequence': getConservationScoresForSubSequence,
        'getTopConservedSpans': getTopConservedSpans,
        'getTopConservedResidues': getTopConservedResidues,
        'compareConservation': compareConservation,
    }

    switcher.get(action, invalid)()

#The main menu defining object
menuAnswer = None
menu = [
    {
        'type': 'list',
        'name': 'choice',
        'message': 'What do you want to do?',
        'choices': [
            {
                'name': 'Load testing files',
                'value': 'loadTestingFiles',
            },
            {
                'name': 'Load FASTA sequence',
                'value': 'loadFasta',
            },
            {
                'name': 'Get the description of a loaded sequence',
                'value': 'getDescription',
            },
            {
                'name': 'Get the sequence of a loaded sequence',
                'value': 'getSequence',
            },
            {
                'name': 'Get the sequence length of a loaded sequence',
                'value': 'getSequenceLength',
            },
            {
                'name': 'Get a subsequence of a loaded sequence',
                'value': 'getSubSequence',
            },
            {
                'name': 'Calculate the hamming distance of two sequences',
                'value': 'getHammingDistance',
            },
            {
                'name': 'Align two sequences using edit distance',
                'value': 'getEditDistanceAlign',
            },
            {
                'name': 'Load structure from PDB file',
                'value': 'loadPdb',
            },
            {
                'name': 'Get structure details',
                'value': 'structureDetails',
            },
            {
                'name': 'Get model details',
                'value': 'modelDetails',
            },
            {
                'name': 'Get chain details',
                'value': 'chainDetails',
            },
            {
                'name': 'Get residue details',
                'value': 'residueDetails',
            },
            {
                'name': 'Get atom details',
                'value': 'atomDetails',
            },
            {
                'name': 'Get atoms within distance of atom',
                'value': 'atomsWithinDistance',
            },
            {
                'name': 'Get residues within distance of atom',
                'value': 'residuesWithinDistance',
            },
            {
                'name': 'Load msa',
                'value': 'loadMsa',
            },
            {
                'name': 'Get number of sequences in MSA',
                'value': 'getMsaSequenceCount',
            },
            {
                'name': 'Get MSA length',
                'value': 'getMsaLength',
            },
            {
                'name': 'Get MSA sequence IDs',
                'value': 'getMsaSequenceIds',
            },
            {
                'name': 'Get MSA sequence',
                'value': 'getMsaSequence',
            },
            {
                'name': 'Get MSA column',
                'value': 'getMsaColumn',
            },
            {
                'name': 'Load scoring matrix',
                'value': 'loadScoringMatrix',
            },
            {
                'name': 'Get MSA column sum of pairs score',
                'value': 'getColumnSumOfPairs',
            },
            {
                'name': 'Get MSA sum of pairs score',
                'value': 'getMsaSumOfPairs',
            },
            {
                'name': 'Get conservation score for sequence in MSA',
                'value': 'getConservationScoreForSequence',
            },
            {
                'name': 'Get sum of pairs score for subsequence of MSA',
                'value': 'getConservationScoresForSubSequence',
            },
            {
                'name': 'Get top conserved spans in MSA',
                'value': 'getTopConservedSpans',
            },
            {
                'name': 'Get top conserved positions in MSA',
                'value': 'getTopConservedResidues',
            },
            {
                'name': 'Compare conservation of active site residues to others',
                'value': 'compareConservation',
            },
            {
                'name': 'Exit',
                'value': 'exit',
            },
        ],
    },
]

#We always clean the screen for the next action. The menu itself is just a loop
if __name__ == '__main__':
    cls()
    while menuAnswer != 'exit':
        menuAnswer = prompt(menu, style=style)['choice']
        doMenuAction(menuAnswer)
        cls()
