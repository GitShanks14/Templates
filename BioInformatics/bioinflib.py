import random

#########################################################################
#  COURSE 0 WEEK ONE FUNCTIONS                                          #
#########################################################################


#Input: A string Text and a string Pattern
#Output: The number of times Pattern appears in Text
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 


#Input: Text and an integer k.
#Output: A dict containing all K-mers present in the string mapped to it's
#        frequency of occurrence.    
def FrequencyMap(Text, k):
    #Initialize the function local variables
    freq = { }
    n = len ( Text )
    
    #Initialize the dict
    for i in range ( n - k + 1 ) :
        Pattern = Text [ i : i + k ]
        freq [ Pattern ] = 0
    
    #Count
    for i in range ( n - k + 1 ) :
        Pattern = Text [ i : i + k ]
        freq [ Pattern ] = freq [ Pattern ] + 1
    
    return freq


# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text

def FrequentWords(Text, k , freq ):
    words = []
    freq = FrequencyMap(Text, k)
    for key in freq:
        if freq[key] >= freq:
            pattern = key
            words.append(pattern)
    return words

# Input : Genome
# Output : Reverse complement of Genome
def ReverseComplement ( Text ) :
    map = { 'A' : 'T' , 'T' : 'A' , 'G' : 'C' , 'C' : 'G' }
    s = ''
    for i in Text [ -1 :  : -1 ] :
        s += map [ i ]
    return s

##########################################################################
# COURSE 0 WEEK TWO                                                      #
##########################################################################


# Input:  Strings Genome and symbol to be searched for in extended genome
# Output: No. of occurences of symbol in every half window of genome
def SymbolArray(Genome, symbol):
    # WARNING: INEFFICIENT VERSION! O(n^2)
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)] , symbol )
    return array


def FasterSymbolArray(Genome, symbol):
    array = { }
    n = len ( Genome )
    ExtendedGenome = Genome + Genome [ 0 : n // 2 ]
    array [ 0 ] = PatternCount ( symbol , Genome [ 0 : n // 2 ] )
    for i in range ( 1 , n ) :
        array [ i ] = array [ i - 1 ]
        if ( ExtendedGenome [ i - 1 ] == symbol ) :
            array [ i ] -= 1
        if ( ExtendedGenome [ i + ( n // 2 ) - 1 ]  == symbol ) :
            array [ i ] += 1   
    return array

# Input : A genome
# Output: The #G encountered - #C encountered when you traverse genome
#           recorded at every index i in the genome. 
def SkewArray(Genome):
    array = [ 0 ]
    n = len ( Genome )
    for i in range ( 0 , n ) :
        if ( Genome [ i : i + 1 ] == 'G' ) : 
            array . append ( array [ i ] + 1 ) 
        elif ( Genome [ i : i + 1 ] == 'C' ) :
            array . append ( array [ i ] - 1 ) 
        else :
            array . append ( array [ i ] )
    return array

# Input: Genome
# Output: minima locations of #G-#C (skew)    
def MinimumSkew(Genome):
    positions = [] # output variable
    skew = SkewArray ( Genome )
    min_val = 0
    for i in range ( len ( Genome ) ) :
        if ( skew [ i ] < min_val ) :
            min_val = skew [ i ]
            positions . clear ( )
            positions . append ( i )
        elif ( skew [ i ] == min_val ) : 
            positions . append ( i )
    return positions


# Input:  Two strings p and q
# Output: Hamming Distance between p and q.
def HammingDistance(p, q):
    h = 0
    a = len ( p ) 
    b = len ( q ) 
    if ( a < b ) :
        l = a
    else : 
        l = b
    for i in range ( l ) :
        if ( p [ i ] != q [ i ] ) :
            h += 1
    return h

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    l = len ( Pattern )
    for i in range ( len ( Text ) - l + 1 ) :
        if ( HammingDistance ( Text [ i : i + l ] , Pattern ) <= d ) :
            positions . append ( i )
    return positions

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    l = len ( Pattern )
    for i in range ( len ( Text ) - l + 1 ) :
        if ( HammingDistance ( Text [ i : i + l ] , Pattern ) <= d ) :
            count += 1
    return count

#########################################################################
# COURSE 0 WEEK 3                                                       #
#########################################################################

# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = { } # initializing the count dictionary
    k  = len ( Motifs [ 0 ] )
    for symbol in "ACGT":
        count[ symbol ] = [ ]
        for j in range ( k ) :
             count [ symbol ] . append ( 0 )
    for i in range ( k ) :
        for j in Motifs :
            count [ j [ i ] ][ i ] += 1
    return count 

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile ( Motifs ) :
    t = len ( Motifs )
    k = len ( Motifs [ 0 ] )
    profile = Count ( Motifs )
    
    for i in "ATGC" :
        for j in range ( k ) :
            profile [ i ][ j ] = profile [ i ][ j ] / t    
    return profile

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus ( Motifs ) :
    k = len ( Motifs [ 0 ] )
    count = Count ( Motifs )
    consensus = ""
    for j in range ( k ) :
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT" :
            if count [ symbol ][ j ] > m :
                m = count [ symbol ][ j ]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return  consensus

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score ( Motifs ) :
    c = Consensus ( Motifs )
    count = Count ( Motifs )
    t = len ( Motifs ) 
    score = 0
    for i in range ( len ( c ) ) : 
        score += t - count [ c [ i ] ][ i ]
    return score

def Pr ( Pattern , profile ) :
    p = 1
    print ( Pattern )
    for i in range ( len ( Pattern ) ) :
        p *= profile [ Pattern [ i ] ][ i ]
    return p


# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
def ProfileMostProbableKmer ( text , k , profile ) :
    p = 0
    s = text [ 0 : k ]
    for i in range ( len ( text ) - k + 1 ) :
        p2 = 1
        for j in range ( k ) :
            p2 *= profile [ text [ i + j ] ][ j ]
        
        if ( p2 > p ) : 
            p = p2
            s = text [ i : i + k ]
    return s


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range ( 0 , t ) :
        BestMotifs . append ( Dna [ i ][ 0 : k ] )
    n = len ( Dna [ 0 ] )
    for i in range ( n - k + 1 ) :
        Motifs = [ ]
        Motifs . append ( Dna [ 0 ][ i : i + k ] )
        for j in range ( 1 , t ) :
            P = Profile ( Motifs [ 0 : j ] )
            Motifs . append ( ProfileMostProbableKmer ( Dna [ j ] , k , P ) )
        if Score ( Motifs ) < Score ( BestMotifs ) :
            BestMotifs = Motifs
    return BestMotifs

##########################################################################
# COURSE 0 WEEK 4                                                        #
##########################################################################

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = { 'A' : [ ] , 'C' : [ ] , 'G' : [ ] , 'T' : [ ] }
    for i in range ( k ) :
        for j in "ACGT" :
            count [ j ] . append ( 1 )
    for i in range ( k ) :
        for j in range ( t ) : 
            count [ Motifs [ j ][ i ] ][ i ] += 1
    return count


# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts ( Motifs )
    for i in range ( k ) :
        for j in "ACGT" : 
            profile [ j ][ i ] /= ( t + 4 ) # +4 as we're adding 1 to A,T,G,C. 
    return profile



#ProfileMostProbableKmerPC ( text , k , profile ) Old version works here too.


def ScorePC ( Motifs ) :
    c = ConsensusPC ( Motifs )
    count = CountWithPseudocounts ( Motifs )
    t = len ( Motifs ) 
    score = 0
    for i in range ( len ( c ) ) : 
        score += t - count [ c [ i ] ][ i ]
    return score

def ConsensusPC ( Motifs ) :
    k = len ( Motifs [ 0 ] )
    count = CountWithPseudocounts ( Motifs )
    consensus = ""
    for j in range ( k ) :
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT" :
            if count [ symbol ][ j ] > m :
                m = count [ symbol ][ j ]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return  consensus

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = [] # output variable
    for i in range ( 0 , t ) :
        BestMotifs . append ( Dna [ i ][ 0 : k ] )
    n = len ( Dna [ 0 ] )
    for i in range ( n - k + 1 ) :
        Motifs = [ ]
        Motifs . append ( Dna [ 0 ][ i : i + k ] )
        for j in range ( 1 , t ) :
            P = ProfileWithPseudocounts ( Motifs [ 0 : j ] )
            Motifs . append ( ProfileMostProbableKmer ( Dna [ j ] , k , P ) )
        if ScorePC ( Motifs ) < ScorePC ( BestMotifs ) :
            BestMotifs = Motifs
    return BestMotifs

# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna). Meant to be used with randomized Profile. 
def Motifs ( Profile , Dna , k ) :
    Motifs = [ ]
    for i in Dna : 
        Motifs . append ( ProfileMostProbableKmer ( i , k , Profile ) )
    return Motifs

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)

# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs ( Dna , k ) :
    pos = [ ]
    for i in range ( len ( Dna ) ) :
        pos . append ( random . randint ( 0 , len ( Dna [ 0 ] ) - k ) )
    Motifs = [ ]
    for i in range ( len ( Dna ) ) :
        Motifs . append ( Dna [ i ][ pos [ i ] : pos [ i ] + k ] ) 
    return Motifs

def RandomizedMotifSearch ( Dna , k , t ) :
    M = RandomMotifs ( Dna , k )
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs ( Profile , Dna , k )
        if Score(M) < Score ( BestMotifs ) :
            BestMotifs = M
        else:
            return BestMotifs 
    
# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    s = 0 
    for i in Probabilities . keys ( ) : 
        s += Probabilities [ i ] 
    for i in Probabilities . keys ( ) : 
        Probabilities [ i ] /= s
    return Probabilities 

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = '' # output variable
    p = random . uniform ( 0 , 1 )
    flag = 1
    for i in Probabilities . keys ( ) :
        if ( p > Probabilities [ i ] ) :
            p -= Probabilities [ i ]
        elif ( flag == 1 ) :
            kmer = i
            flag = 0
    return kmer

#Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len ( Text )
    probabilities = { } 
    for i in range ( 0 , n - k + 1 ) :
        probabilities [ Text [ i : i + k ] ] = Pr ( Text [ i : i + k ] , profile )
    probabilities = Normalize ( probabilities )
    return WeightedDie ( probabilities )

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler ( Dna , k , t , N ) :
    BestMotifs = [ ] # output variable
    M = RandomMotifs ( Dna , k )
    BestMotifs = M
    profile = ProfileWithPseudocounts ( M )
    for i in range ( N ) :
        j = random . randint ( 0 , t - 1 )
        text = ProfileGeneratedString ( M [ j ] , profile , k )
        M [ j ] = ProfileGeneratedString ( text , profile , k )
        profile = ProfileWithPseudocounts ( M )
        if Score ( M ) < Score ( BestMotifs ) :
            BestMotifs = M
    return BestMotifs


##########################################################################
# TESTING GROUNDS                                                        #
##########################################################################


def PatternPositions ( Text , Pattern ) :
    positions = [ ] 
    k = len ( Pattern )
    for i in range ( len ( Text ) - k + 1 ) :
        if ( Pattern == Text [ i : i + k ] ) : 
            positions . append ( i )
    return positions

Text = 'CAAAATTGTCACAACGCTGTGTTGCTTAAGACGATATTTCGGGTAACATGTTTGCGGCCAAGTTAGTGCGCAGCGGCGTTGAGTTCCATCAGGAAAGGATGCTGGTGCCCCGACACGTTGGCTCCCCTGGCGTTAATAAGCTTATAAATGCCGTGGACCGTGTATTAATAAACAGTTCCGCGCGTACAGCCATGAAAAATGCAAGGAGTGCGATAACAATAAACAGAAGTTGTTAGTCTAGGAAGACTATTTTGAGCAGATACGCCCAAGTCTAGGTGTATAAAATCCGAGCTAGATGTCATCCATGGACTATCATTTTACAGTAATTCAAACCCAAATACCCCATACAGGTGACGCGACATCTCCTTTAAGACAGACTTTAAGACTTGACATATACATACACACATAATGCGTTCCGGCTGACTCCTCTCCATGTGCCATGTGGCCCTACACCGTCCTGGTTAGATTTAAAGCATTTGCCAACCCCACACCACACCCACACTCGGACACTCGGCTAACCGCTGAACAGGTCGTTCCCGATGTGGGATGCAACTGGCATAATAAACTAACGGGTCAGGGATTCTACGAGAGTGATCCCTTCGAGAAGATCAATGGAGCCACTTGTAAGCACTGTATAGACGTTCGCAATTGATCGCAATTGACGCAATTGAATGGTCTCCGACTGCTCCCATGTATTTCGACCCCAATGCGGCGTTCCTACCACGTACACACCACGTGAGAGAGCCAGGCCATCCAAATCAGGCACTTGGCCTAATTGGTTGAGCGGTTAGCCGAACCGGACGTGCCGCGCACCGCGCACGCGCACCAACACCAAAACTTAGCGTAGGGCGTTTATGGTGCAAGGAAACTACGCCTCGAAGGGTTGAATGTCCGTTCGTCCGAGAGAACCACAGGGGCGGCTCTAAATCGGCTCTAAATAAAACAAAACCGATCAACCTCGGTCAGTGCTTCCCAAACAACAGAGTCAGAGTGCAGAGTGACAGAGTGAGACTGCGAACTAGACGAACTAGAAGACTAGAAGGGTTGAACTCGTGGCACGCGCCCGGTCCAGCCCGAGTATCTAGCTACGATACCCTCTGTGGAAGTCTGTCTGATCGCTCCAGGTGACGCTCAAACACGTTATCGGTCGTGTAAATTCATCCCCCAACAACAATTTATAATTTATTCCAATTTATTCAGGCTATGATGCCGGGCATAATTGCCAGTATAATCCTGTTTGCTATTTTATAGGATGGAGCCGATTGAATGGCGATCGTTCGAATCGTCGTTACTGAACAGCAATAAAAATTGCCGTCAATTAGGCTGACTTTCCTCAGGGCGAGCTTATGTAACGTCGCTGAAGGATAATTTAGTAAGAGTCGTCGCACCACGTAAGCTCAGTGGTTAGCAGTGAAATAAATGAAATAAAGGTGAAATAACTAATTACCTGGGCGTTCGTTTCCGGCGACATCGGGCGCTTATGCGCCACGACAGCCAATCTATGCGCACAAGGTCAAACCTAGTCTGCTCGCCGGCGCAATTCACTTTCAATACCTTTCAATGTTGGAAGTTATGTATGTTGGAATGAGCTCGCACAAGGTCGATCCACCGGAGCAAAAGAAACGCGAAAGCTTTAAGTTCAGATAAACGCGAGGGCCGCAATGTACTGCATGTTATGTGGTCAATTTAATTTAATTTCGGGTGCCACGGACGCACCCTAGTGAAGCGCCGAGATCTTGGCTCGCGATTTTAGATTAGTTAGATGCATGTTAGATGCAATATCGCACTGGGATGCCTGAATATCCTGTCACGACCGCGAGTCAAGCCATGCCAGGACGCCGCCCCTAATGGATCATAAAATGGGAACAATGGGTCGGGGGTCTCGGGGGTCGGGTCGGTCTCCAATGGGAACAATGGGAACAATGGGAACCTCGGTACGCTCGGTACGCTCGGTACGCTCCTCCCTCGGTACGTACG'
L = 29
k = 9
freq = 4  

w = FrequentWords(Text, k , freq )
for i in w :
    print ( i , end = ' ' )



















