def hammingDistance(sequenceA, sequenceB):
    '''
    Computes the hammingDistance between two sequences

            Parameters:
                sequenceA (Sequence): The first sequence to align
                sequenceB (Sequence): The second sequence to align

            Returns:
                hammingDistance (int): The hamming distance
    '''
    if(sequenceA.getSequenceLength() != sequenceB.getSequenceLength()):
        raise ValueError('Sequences must be of the same length to be comparable using hamming distance!')

    #Hamming distance is calculated by iterating over the sequences and counting mismatches
    differences = 0
    for i in range(0, sequenceA.getSequenceLength()):
        if sequenceA.getSequence()[i] != sequenceB.getSequence()[i]:
            differences += 1

    return differences
