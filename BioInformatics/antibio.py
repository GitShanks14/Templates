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
   
