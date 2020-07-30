# Bioinformatics toolbox

A toolbox created for a MFF.CUNI.CZ course.

## Description

Implementations of basic methods from the field of sequence and structural bioinformatics such as parsing files or assessing similarity of sequences and structures.

### Namely the tasks:

_The tasks presume you have completed the required loading steps or used the Shortcut loading for testing purposes as described in its section, otherwise an error message is shown_

* [Processing FASTA files](http://bioinformatika.mff.cuni.cz/repository/#/assignments/detail?id=tools_fasta_parsing)
  1. Read in a FASTA file with an arbitrary number of molecules
    * Example menu path
      1. Select **Load FASTA sequence**
      1. Enter **TestingFiles/A9WZ33Double.fasta**
      1. Enter an integer specifying which sequence from the file you want
  1. Obtain description/sequence of any of the molecules.
    * Example menu path
      1. Select **Get the description of a loaded sequence** or **Get the sequence of a loaded sequence**
      1. Select the sequence you wish to display the data for
  1. Return sequence length for given sequence.
    * Example menu path
      1. Select **Get the sequence length of a loaded sequence**
      1. Select the sequence you wish to display the data for
  1. Return subsequence of given sequence.
    * Example menu path
      1. Select **Get a subsequence of a loaded sequence**
      1. Select the sequence you wish to display the data for
      1. Enter the beginning index integer
      1. Enter the end index integer
* [Measuring sequence similarity using Hamming distance](http://bioinformatika.mff.cuni.cz/repository/#/assignments/detail?id=tools_hamming_distance)
  * Example menu path
    1. Select **Calculate the hamming distance of two sequences**
    1. Select the first sequence
    1. Select the second sequence  
* [Sequence alignment using edit distance](http://bioinformatika.mff.cuni.cz/repository/#/assignments/detail?id=tools_edit_distance)
  * Example menu path
    1. Select **Align two sequences using edit distance**
    1. Select the first sequence
    1. Select the second sequence  
* [Processing PDB files](http://bioinformatika.mff.cuni.cz/repository/#/assignments/detail?id=tools_pdb_parsing)
  1. Read in a PDB file.
    * Example menu path
      1. Select **Load structure from PDB file**
      1. Enter **TestingFiles/3owe.pdb**
  1. Obtain an object representing a model
    * Example menu path
      1. Select **Get model details**
      1. Select the structure
      1. Select the model
  1. Obtain an object representing a structure (chain) within a model.
    * Example menu path
      1. Select **Get chain details**
      1. Select the structure
      1. Select the model
      1. Select the chain
  1. Obtain an object representing residuum within a chain.
    * Example menu path
      1. Select **Get residue details**
      1. Select the structure
      1. Select the model
      1. Select the chain
      1. Select the residue
  1. Obtain an object representing an atom within a residue
    * Example menu path
      1. Select **Get atom details**
      1. Select the structure
      1. Select the model
      1. Select the chain
      1. Select the residue
      1. Select the atom
  1. Obtain information about the stored structure (number of models, structures, residues, atoms).
    * Example menu path
      1. Select **Get structure details**
      1. Select the structure      
  1. Compute the width of the structure (maximum of distance of any two atoms).
    * Is a part of the above mentioned choice **Get structure details**   
  1. Obtain list of atoms being in given distance from given ligand (HETATM).
    * Example menu path  
      1. Select **Get atoms within distance of atom**
      1. Select the structure
      1. Select the model
      1. Select the chain
      1. Select the residue
      1. Select the atom
      1. Enter the distance
  1. Obtain list of residues being in given distance from given ligand (HETATM).
    * Example menu path   
      1. Select **Get residues within distance of atom**
      1. Select the structure
      1. Select the model
      1. Select the chain
      1. Select the residue
      1. Select the atom
      1. Enter the distance    
* [Processing multiple sequence alignment](http://bioinformatika.mff.cuni.cz/repository/#/assignments/detail?id=tools_msa_parsing)
  1. Read and parse MSA.
    * Example menu path  
      1. Select **Load msa**
      1. Enter **TestingFiles/p53.clustal**
  1. Retrieve sequence by its position or ID.
    * Example menu path   
      1. Select **Get MSA sequence**
      1. Select the MSA
      1. Select the sequence
  1. Retrieve given column from the MSA
    * Example menu path   
      1. Select **Get MSA column**
      1. Select the MSA
      1. Enter the column number
  1. Retrieve sum of pairs score of a column and whole MSA with respect to given scoring matrix.
    * Example menu path   
      1. Select **Load scoring matrix**
      1. Enter **TestingFiles/Blosum62**
      1. Select **Get MSA column sum of pairs score** or **Get MSA sum of pairs score**
      1. Select the MSA
      1. Select the scoring matrix
      1. If you selected the column variant, enter the column number
* [Conservation determination from multiple aligned sequences](http://bioinformatika.mff.cuni.cz/repository/#/assignments/detail?id=tools_msa_properties)
  * Example menu path
    1. Select **Get top conserved spans in MSA** or **Get top conserved positions in MSA**
    1. Select the MSA
    1. Select the scoring matrix
    1. Enter the number of spans or positions you want to get
* [Combining structure and sequence](http://bioinformatika.mff.cuni.cz/repository/#/assignments/detail?id=tools_structure_sequence)
  * Example menu path
    1. Select **Compare conservation of active site residues to others**
    1. Select structure
    1. Select MSA
    1. Select scoring matrix
    1. Select chain
    1. Enter atom type to search around
    1. Enter the distance you want to search in
  * Results:
    * Up to distance 2 we have only HOH and the ZN itself, then we get more interesting results. With distance 3 we get the 3 CYS that are a part of the active site
    * Why is position 5 greater? Because it is also conserved and contains a higher scored residue (W-W score is 11 while C-C is 9 in our Blosum testing SM)
    * With the above point explained, the results were as expected - Active site was more conserved than the other parts of the protein

## Instalation

* Please install Python v 3.7.4 and pip package manager
* In the directory containing this program create a python virtual environment
  ~~~
  python -m venv env
  ~~~
* Activate said virtual environment
  * On windows:
    ~~~
    .\env\Scripts\activate
    ~~~
  * On unix:
    ~~~
    source env/bin/activate
    ~~~
* Install required packages using the command:
  ~~~
  pip install -r requirements.txt
  ~~~
  * Notable packages include:
    * PyInquirer
      * Menu handling
    * biopython
      * Handling of bioinformatic files
    * sortedcontainers
      * Sorted container implementations
    * scipy
      * Used for 3D convex hull computation

## Directory structure

* docs
  * Contains pydoc documentation for the toolbox
* Exceptions
  * Contains custom exception classes
* SequenceTools
  * Contains classes that deal with sequences
* StructureTools
  * Contains classes that deal with structures
* menu.py
  * The base menu used for running the package from command line

## Start

Write the following into command line when you are in the directory containing this package and have the installed environment active.

~~~
python menu.py
~~~

## Stop

Select **Exit** in menu

## Shortcut loading for testing purposes

Choose the menu field: **Load testing files** then enter **1** or **2** to select sequence as one of the testing fasta files contains multiple sequences
