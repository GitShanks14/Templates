###############################################################################
#                        Motif finder                                         #
###############################################################################

import numpy as np

######################### Utility functions ###################################

# Input : Genome
# Output : Reverse complement of Genome
def ReverseComplement ( Text ) :
    m = { 'A' : 'T' , 'T' : 'A' , 'G' : 'C' , 'C' : 'G' }
    s = ''
    for i in Text [ -1 :  : -1 ] :
        s += m [ i ]
    return s

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

# Input  : A pattern and a distance d
# Output : All patterns in the d neighborhood 

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

# Brute force Motif search
# Input  : collection of Dna strings , k , d
# Output : All motifs of len k that have at least 1 kmer from their d neighbors
            # in each string of Dna
        
def MotifSearchBF ( dna , k , d ) :
    patterns = [ ]
    for i in range ( 0 , len ( dna [ 0 ] ) - k + 1 ) :
        neighbors = Neighbors ( dna [ 0 ][ i : i + k ] , d )
        for j in neighbors:
            count = 0
            for l in dna :
                for i in range ( 0 , len ( l ) - k + 1 ) :
                    if HammingDistance ( j , l [ i : i + k ] ) <= d :
                        count += 1
                        break
            if count == len ( dna ) :
                patterns . append ( j )
    Patterns = [] 
    [ Patterns . append ( x ) for x in patterns if x not in Patterns ] 
    return Patterns

# Input:  A set of kmers Motifs
# Output: Counts the no. of each bases in each position naively w/o pseudocounts
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

# Input:  A set of kmers Motifs
# Output: Counts the no. of each base in each position, with a pseudocount of 1
def CountPC ( Motifs ) :
    t = len ( Motifs )
    k = len ( Motifs [ 0 ] )
    count = { 'A' : [ ] , 'C' : [ ] , 'G' : [ ] , 'T' : [ ] }
    for i in range ( k ) :
        for j in "ACGT" :
            count [ j ] . append ( 1 )
    for i in range ( k ) :
        for j in range ( t ) : 
            count [ Motifs [ j ][ i ] ][ i ] += 1
    return count

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus ( count ) :
    k = len ( count [ 'A' ] )
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
# Output: The naive score of these k-mers without pseudocount
def Score ( count ) :
    c = Consensus ( count )
    t = 0
    for i in "ACGT" :
        t += count [ i ][ 0 ]
    score = 0
    for i in range ( len ( c ) ) : 
        score += t - count [ c [ i ] ][ i ]
    return score

# ScoreV2 is WORSE than V1. timed to be 20% SLOWER than Score. AVOID. 
def ScoreV2 ( Motifs ) :
    c = Consensus ( Count ( Motifs ) )
    score = 0
    for i in range ( len ( Motifs ) ) :
        score += HammingDistance( Motifs [ i ] , c )
    return score

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile ( Motifs ) :
    t = len ( Motifs )
    k = len ( Motifs [ 0 ] )
    count = Count ( Motifs )
    profile = np . zeros ( ( 4 , k ) )
    m = { 'A' : 0 , 'C' : 1 , 'G' : 2 , 'T' : 3 }
    for i in "ATGC" :
            profile [ m [ i ] ] = np . array ( count [ i ] ) / t    
    return profile    

def ProfilePC ( Motifs ):
    t = len ( Motifs )
    k = len ( Motifs [ 0 ] )
    count = CountPC ( Motifs )
    profile = np . zeros ( ( 4 , k ) )
    m = { 'A' : 0 , 'C' : 1 , 'G' : 2 , 'T' : 3 }
    for i in 'ATGC' :
            profile [ m [ i ] ] = np . array ( count [ i ] ) / ( t + 4 ) 
    return profile
    
def Entropy ( profile ) :
    p = - profile * np . log2 ( profile + 0.00000000001 )
    p = np . sum ( p )
    return p

# Input  : A pattern and a list of DNA strings of same length 
# Output : The best motif set built around the pattern is selected and returned
#           The score and Motif set are returned as a tuple. 
# O ( n.t.k )
def BestMotifs ( Pattern , Dna ) :
    Motifs = [ ]
    score = 0
    k = len ( Pattern )
    for i in range ( len ( Dna ) ) :                    # t steps
        d = 100000
        p = None
        for j in range ( len ( Dna [ 0 ] ) - k + 1 ) :  # n - k ~ n steps
            d2 = HammingDistance ( Dna [ i ][ j : j + k ] , Pattern ) # k steps
            if ( d2 < d ) :
                d = d2
                p = Dna [ i ][ j : j + k ]
        Motifs . append ( p )
        score += d 
    return score , Motifs
            
# BRUTE FORCE MOTIF SEARCH: Traverses ALL possible kmers to find the one for 
# which the lowest scoring motif set exists. DOUBLE MINIMISATION PROBLEM
# Input : A list of strings of Dna of equal length, length k of pattern to be found
# INEFFICIENT. 4^k can be narrowed down dramatically. 
# Method may not be improved upon, due to existence of more promising alternative.
def BFMotifSearch ( Dna , k ) :
    Motifs = [ ]
    p = [ ]
    d = 10000000
    for i in range ( 4 ** k ) :
        p2 = NumberToPattern( i , k )
        d2 , m2 = BestMotifs ( p2 , Dna )
        if ( d2 < d ) :
            p . clear ( )
            Motifs . clear ( )
            p . append ( p2 )
            Motifs . append ( m2 )
            d = d2
        elif ( d2 == d ) :
            p . append ( p2 )
            Motifs . append ( m2 )
    return d , Motifs , p
            
###############################################################################
#                            TESTING GROUNDS                                  #
###############################################################################

import time
t0 = time . time ( )

def timer ( ):
    global t0
    print ( "Execution of block took {} s" . format ( time . time ( ) - t0 ) )
    t0 = time . time ( )
    

k = 6

with open ( '/Users/sashank/Desktop/Data/d.txt' , mode = 'r') as f:
    genome = f . readline ( 10000000000000000 )
    Dna = genome . split ( )
    print ( "Data read from file." )
    score , motifs , p = BFMotifSearch ( Dna , k )
    print ( score ) 
    for i in p :
        print ( i , end = ' ' )
    print ( '' )
    timer ( )
    
#### DOCUMENT THE CODE!!!

    