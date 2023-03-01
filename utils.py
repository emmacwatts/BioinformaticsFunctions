def plaintextFASTA(header, text, filepath = 'auto', nucleotide = True, allPotentialGenes = False):
    """ Make a formatted FASTA file from unformatted bases (e.g. includes spaces, 
    position numbers, etc.). Will also check for lack of start/stop, non-codon-divisible, false bases/amino acids.
    If looking for just coding sequence, can mount to start at ATG (atgStart = True).
    

    params:
        filepath: file path and name for output file
        header: desired header for the fasta file
        text: the sequence to be formatted
        nucleotide: If true, the sequence will be treated as a nucleotide sequence. If false, will be treated as a protein sequence.
        allPotentialGenes: If true, will find all potential genes from the text (for all frames, any ATG until stop).
    """
    #Define the automatic filepath (the header becomes the filename)
    if filepath == 'auto':
        filepath = f"{header}.fasta"

    #Define the text type
    if nucleotide == True:
        valid = ['A', 'T', 'G', 'C', 'N'] #N for any base
    else:
        valid = ['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    #Ensure str is uppercase
    text = text.upper()

    #Error messages start as false
    removedNumbers = False
    removedInvalidCharacters = False
    removedSpecialCharacter = False
    invalidCharacters = []

    #Format the text - removes any non-valid characters, numbers, etc.
    formattedText = ""
    for character in text:
        if character in valid:
            formattedText += character
        if character.isnumeric() == True:
            removedNumbers = True
        if character.isalpha() and character not in valid:
            removedInvalidCharacters = True
            invalidCharacters.append(character)
        if character.isnumeric() == False and character.isalpha() == False:
            removedSpecialCharacter = True

    #Error messages - warn about what has been removed from original text
    if len(formattedText) %3 != 0:
        print("Invalid codons: number of bases not divisible by 3.")
    if removedNumbers == True:
        print("Removed numbers from original sequence.")
    if removedInvalidCharacters == True:
        print(f"Removed invalid characters {invalidCharacters}.")
    if removedSpecialCharacter == True:
        print("Removed special charactes.")


    #Write fasta file
    with open(filepath, 'a') as f:
        f.write(f"\n>{header}\n{formattedText}")

def baseContent(dnaFASTA, GCcontent = False):
    """
    For a DNA fasta file, calculates content of A, T, G, C. 
    """

    from Bio import SeqIO

    baseCounts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N':0}

    for fasta in SeqIO.parse(dnaFASTA, "fasta"):
        currentProtein = fasta.id
        sequence = fasta.seq
        for base in sequence:
            if base not in baseCounts.keys():
                print(f"Non-standard base in {currentProtein}, '{base}'.")
            else:
                baseCounts[base] += 1

    if GCcontent == True:
        GC = (baseCounts['G'] + baseCounts['C'])/sum(baseCounts.values())
        return GC
    else:
        return baseCounts
    


