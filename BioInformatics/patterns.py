###############################################################################
#  Finding Patterns                                                           #
###############################################################################

# Input : Genome
# Output : Reverse complement of Genome
def ReverseComplement ( Text ) :
    m = { 'A' : 'T' , 'T' : 'A' , 'G' : 'C' , 'C' : 'G' }
    s = ''
    for i in Text [ -1 :  : -1 ] :
        s += m [ i ]
    return s

def PatternMatching ( Pattern , Genome ) :
    pos = [ ]
    l = len ( Pattern )
    for i in range ( len ( Genome ) - l + 1 ) :
        if Genome [ i : i + l ] == Pattern :
            pos . append ( i )
    return pos

def PatternToNumber ( text ) :
    s = 0
    m = { 'A' : 0 , 'C' : 1 , 'G' : 2 , 'T' : 3 } 
    for i in text :
        s *= 4
        s += m [ i ]
    return s

def NumberToPattern(index, k):
    m = { 0 : 'A' , 1 : 'C' , 2 : 'G' , 3 : 'T' }
    s = ''
    while k :
        k -= 1
        s += m [ index % 4 ]
        index //= 4
    return s [ -1 : : -1 ]

def ComputingFrequencies(text, k):
    freq = [ ]
    for i in range ( 4 ** k ) :
        freq . append ( 0 )
    for i in range ( len ( text ) - k + 1 ) :
        freq [ PatternToNumber ( text [ i : i + k ] ) ] += 1
    return freq

################################################################################

#Input: A string Text and a string Pattern
#Output: The number of times Pattern appears in Text
def PatternCount ( Text , Pattern ) :
    count = 0
    for i in range ( len ( Text ) - len ( Pattern ) + 1 ):
        if ( Text [ i : i + len ( Pattern ) ] == Pattern ) :
            count += 1
    return count 


#Input: Text and an integer k.
#Output: A dict containing all K-mers present in the string mapped to it's
#        frequency of occurrence.    
def FrequencyMap(Text, k):
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


# Input  :  A string Text and an integer k
# Output : A list containing all k-mers that appear more than f times in Text.
def FrequentWords ( Text , k , f ) :
    words = [ ]
    freq = FrequencyMap ( Text , k )
    for key in freq :
        if freq [ key ] >= f :
            pattern = key
            words . append ( pattern )
    return words

# Input  : Input: A string Genome, and integers k, L, and t.
# Output : All distinct k-mers forming ( L , t ) clumps in Genome.
# An ( L , t ) clump is a clump that appears at least t times in a window of
#   Length L.

def LtClumps ( Text , L , t , k ) :
    clumps = set ( [ ] )
    for i in range ( len ( Text ) - L + 1 ) :
        words = FrequentWords ( Text , k , t )
        clumps . update ( words )
    return clumps

#This approach is inefficient.
#Better approach: make dict mapping every kmer in Text to a list containing
#                 the indices where the kmer occurs. O( |Text| )
#                 Now, scan through the list corresponding to each kmer and 
#                 If there are t indices that fall inside a box of size L-k,
#                 then the said kmer is an L,t-Clump. 


# Input  : Genome and kmer length.
# Output : a dict mapping kmers to a list containing the indices of occurences

def KmersToPositions ( Text , k ) :           # NET: O ( |Text| )
    m = { }
    for i in range ( len ( Text ) - k + 1 )  :  # O ( |Text| )
        if ( Text [ i : i + k ] in m ) :
            m [ Text [ i : i + k ] ] . append ( i )
        else :
            m [ Text [ i : i + k ] ] = [ i ]
    return m

# Input  : Genome, Window length, frequency cutoff, and kmer length.
# Output : List of Kmers that are (L,t)-clumps. 
# O ( |Text| ) as that's the upper bound on no. of indices.
def LtClumpsV2 ( Text , L , t , k ) : 
    m = KmersToPositions ( Text , k )
    clumps = [ ]
    for key, l in m . items ( ) :           
        if ( len ( l ) - t + 1 <= 0 ) :
            continue
        for i in range ( len ( l ) - t + 1 ) : 
            if ( l [ i + t - 1 ] - l [ i ] + k <= L ) :
                clumps . append ( key )
                break
    return clumps


###############################################################################
#                       FINDING THE ORI                                       #
###############################################################################

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
def ApproximatePatternCount(Text, Pattern, d):
    count = 0 # initialize count variable
    l = len ( Pattern )
    for i in range ( len ( Text ) - l + 1 ) :
        if ( HammingDistance ( Text [ i : i + l ] , Pattern ) <= d ) :
            count += 1
    return count

def ApproximateFrequentWords ( Text , k , d ) :
    words = [ ]
    freq  = { }
    m = 0
    # Creating dict of all possible kmers
    for i in range ( 4 ** k ) :
        freq [ NumberToPattern ( i , k ) ] = 0
    #Now, for each key in freq, count the number of approximate matches and
    # simultaneously find the max
    for key , value in freq.items ( ) :
        value = ApproximatePatternCount ( Text , key , d )
        if ( value > m ) :
            m = value
            words . clear ( )
            words . append ( key )
        elif ( value == m ) :
            words . append ( key )
    return words

def ApproximateFrequentWordsRC ( Text , k , d ) :
    words = [ ]
    freq  = { }
    m = 0
    
    # Creating dict of all possible kmers
    for i in range ( 4 ** k ) :
        s = NumberToPattern ( i , k )
        freq [ s ] = ApproximatePatternCount ( Text , s , d )
    print ( "Dict Created" , time . time ( ) - t1 )
    # Now total the scores and find max
    for key in freq :
        v = freq [ key ] + freq [ ReverseComplement ( key ) ]
        if ( v > m ) :
            m = v
            words . clear ( )
            words . append ( key )
        elif ( v == m ) :
            words . append ( key )
    return words

def Neighbors(Pattern, d):
    if d == 0 :
        return [ Pattern ]
    if len ( Pattern ) == 1 : 
            return [ 'A' , 'C' , 'G' , 'T' ]
    Neighborhood = [ ]
    bases = [ 'A' , 'C' , 'G' , 'T' ]
    SuffixNeighbors = Neighbors ( Pattern [ 1 : ] , d )
    for Text in SuffixNeighbors :
        if HammingDistance ( Pattern [ 1 : ] , Text ) < d :
            for x in bases :
                    Neighborhood . append ( x + Text )
        else :
            Neighborhood . append ( Pattern [ 0 ] + Text )
    return Neighborhood



    
    
    
    
    
    
    
    
    
