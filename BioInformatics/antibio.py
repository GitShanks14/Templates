######################### SEQUENCING ANTIBIOTICS ##############################

# Input  : RNA
# Output : Protein translated to, if read assuming 0th reading frame
def DnaToRna ( Dna ) :
    Rna = ''
    for i in Dna :
        if ( i != 'T' ) :
            Rna += i
        else :
            Rna += 'U' 
    return Rna

# Input : Genome
# Output : Reverse complement of Genome
def ReverseComplement ( Text ) :
    m = { 'A' : 'T' , 'T' : 'A' , 'G' : 'C' , 'C' : 'G' }
    s = ''
    for i in Text [ -1 :  : -1 ] :
        s += m [ i ]
    return s

def ReverseComplementRna ( Text ) :
    m = { 'A' : 'U' , 'U' : 'A' , 'G' : 'C' , 'C' : 'G' }
    s = ''
    for i in Text [ -1 :  : -1 ] :
        s += m [ i ]
    return s

def RnaTranslation ( Rna ) :
    m = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    n = len ( Rna )
    Protein = ''
    i = 0
    while ( i < n - 2 ) :
        AminoAcid = m [ Rna[ i : i + 3 ] ]
        if AminoAcid != 'STOP' : 
            Protein += AminoAcid
        else :
            return Protein
        i += 3
    return Protein

# Input  : Dna , Peptide
# Output : List of all substrings that map to the Peptide
# O ( nk )
def FindAllEncodings ( Dna , Peptide ) :
    Rna = DnaToRna( Dna )
    k = len ( Peptide ) * 3
    n = len ( Dna ) 
    substrings = [ ]
    for i in range ( n - k + 1 ) : 
        s = Rna [ i : i + k ]
        RNA = s
        RNARC = ReverseComplementRna ( s )
        if RnaTranslation ( RNA ) == Peptide or RnaTranslation ( RNARC ) == Peptide :
            substrings . append ( Dna [ i : i + k ] ) 
    return substrings
   
mass = {
        'G' : 57,
        'A' : 71,
        'S' : 87,
        'P' : 97,
        'V' : 99,
        'T' : 101,
        'C' : 103,
        'I' : 113,
        'L' : 113,
        'N' : 114,
        'D' : 115,
        'K' : 128,
        'Q' : 128,
        'E' : 129,
        'M' : 131,
        'H' : 137,
        'F' : 147,
        'R' : 156,
        'Y' : 163,
        'W' : 186
}

def LinearSpectrum ( Peptide ) :
    PrefixMass = [ 0 ]
    Spectrum = [ 0 ]
    n = len ( Peptide )
    for i in range ( n ) :
        PrefixMass . append ( PrefixMass [ i ] + mass [ Peptide [ i ] ] )
    for i in range ( n ) :
        for j in range ( i + 1 , n + 1 ) :
            Spectrum . append ( PrefixMass [ j ] - PrefixMass [ i ] )
    Spectrum . sort ( )
    return Spectrum

def CyclicSpectrum ( Peptide ) :
    PrefixMass = [ 0 ]
    Spectrum = [ 0 ]
    n = len ( Peptide )
    Peptide += Peptide
    for i in range ( 2 * n ) :
        PrefixMass . append ( PrefixMass [ i ] + mass [ Peptide [ i ] ] )
    for i in range ( n ) :
        for j in range ( i + 1 , i + n ) :
            Spectrum . append ( PrefixMass [ j ] - PrefixMass[ i ] )
    Spectrum . sort ( )
    Spectrum . append ( PrefixMass [ n ] )
    return Spectrum


CountMemoize = { 0 : 1 }
AminoAcids = [ 'G' , 'A' , 'S' , 'P' , 'V' , 'T' , 'C' , 'I' , 'N' , 
               'D' , 'K' , 'E' , 'M' , 'H' , 'F' , 'R' , 'Y' , 'W' ]

def CountPeptidesOfMass ( m ) : #Does NOT account for amino acids of same mass
    count = 0
    if ( m in CountMemoize ) :
        return CountMemoize [ m ]
    elif ( m < 0 ) :
        return 0
    for i in AminoAcids :
        x = m - mass [ i ]
        c = CountPeptidesOfMass ( x )
        count += c
        if x not in CountMemoize :
            CountMemoize [ x ] = c
    return count
        
def CycloPeptideSequencing ( Spectrum ) :
    Candidates = [ '' ]
    pmass = Spectrum [ -1 ]
    FinalPeptides = [ ]    
    while ( Candidates ) :
        # Grow
        NewCandidates = [ ]
        for i in Candidates :
            for j in AminoAcids :
                NewCandidates . append ( str ( i ) + str ( j ) )
        Candidates = NewCandidates . copy ( )
        print ( len ( Candidates ) , ' ->' , end = ' ' )
        # Filter
        for i in NewCandidates : 
            sp = LinearSpectrum ( i )
            notin = False
            for j in sp :
                if j not in Spectrum :
                    Candidates . remove ( i )
                    notin = True
                    break
            if ( notin == False and sp [ -1 ] == pmass ) :
                if ( i not in FinalPeptides ) :
                    FinalPeptides . append ( i )
                    Candidates . remove ( i )
        print ( len ( Candidates ) ) 
        timer ( )
    return FinalPeptides 

