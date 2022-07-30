# Author: <Davin Marro> <dmarro@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
Bi621 Additions up to v.0.6 completed 7/17/22'''

__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

from asyncio import base_subprocess
from email.mime import base


DNA_bases = "ATCG"
RNA_bases = "ATCGU"

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    qual_score = (ord(letter) -33)
    return qual_score

def qual_score(phred_score: str) -> float:
    score_list = []
    """Averages the quality score of a string using convert_phred"""    
    for phred in phred_score:
        score = convert_phred(phred)
        score_list.append(score)
    average = sum(score_list) / len(score_list)
    return average #message = "The average quality score of this string is:"

def validate_base_seq(base_sequence, RNAflag = False):
    '''Determines whether string(seq) is a DNA sequence (ATCG)'''
    seq = base_sequence.upper()
    return len(seq) == (seq.count('U' if RNAflag else 'T') +
                        seq.count('C') +
                        seq.count('A') +
                        seq.count('G'))
def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

def oneline_fasta():
    '''Compiles fasta sequence lines into one line'''
    with open('input.fasta') as f_input, open('output.fasta', 'w') as f_output:
        block = []
    for line in f_input:
        if line.startswith('>'): #Targets header lines
            if block: 
                f_output.write(''.join(block) + '\n') #Joins sequences inclusive of newline char.
                block = []
            f_output.write(line)
        else:
            block.append(line.strip())

    if block:
        f_output.write(''.join(block) + '\n')
    pass


def validate_DNA_seq(DNA):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts, Gs, and Cs. False otherwise. Case insensitive.'''
    #here is a comment
    DNA = DNA.upper()
    return len(DNA) == DNA.count("A") + DNA.count("T") + DNA.count("G") + DNA.count("C")

if __name__ == "test.fa":
    #Unit test fasta for functions
    pass
