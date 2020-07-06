###############################################################################
#                        Motif finder                                         #
###############################################################################

import numpy as np
import random 
import sys

# Edit the path so as to import BioUtils
sys . path . insert ( 0 , '/Users/sashank/Desktop/Python/Libraries/Bioinformatics/BioUtils' )

import BioUtils as bf



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
        score += bf . HammingDistance( Motifs [ i ] , c )
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
            d2 = bf . HammingDistance ( Dna [ i ][ j : j + k ] , Pattern ) # k steps
            if ( d2 < d ) :
                d = d2
                p = Dna [ i ][ j : j + k ]
        Motifs . append ( p )
        score += d 
    return score , Motifs
            
def Pr ( string , profile ) :
    m = { 'A' : 0 , 'C' : 1 , 'G' : 2 , 'T' : 3 }
    p = 1
    for i in range ( len ( string ) ) :
        p *= profile [ m [ string [ i ] ] , i ]
    return p

def MostProbableKmer ( genome , profile , k ) :
    p = 0
    kmer = [ ]
    for i in range ( len ( genome ) - k + 1 ) :
        p2 = Pr ( genome [ i : i + k ] , profile )
        if ( p2 > p ) : 
            p = p2
            kmer . clear ( )
            kmer . append ( genome [ i : i + k ] )
        elif ( p2 == p ) :
            kmer . append ( genome [ i : i + k ] )
    return kmer [ 0 ]

# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna). Meant to be used with randomized Profile. 
def Motifs ( Profile , Dna , k ) :
    Motifs = [ ]
    for i in Dna : 
        Motifs . append ( MostProbableKmer ( i , Profile , k ) )
    return Motifs

def RandomMotifs ( Dna , k ) :
    t = len ( Dna )
    n = len ( Dna [ 0 ] )
    Motifs = [ ]
    for i in range ( t ) :
        x = random . randint ( 0 , n - k )
        Motifs . append ( Dna [ i ][ x : x + k ] ) 
    return Motifs

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


###############################################################################
#                       MOTIF SEARCHES                                        #
###############################################################################

# Brute force Motif search
# Input  : collection of Dna strings , k , d
# Output : All motifs of len k that have at least 1 kmer from their d neighbors
            # in each string of Dna
        
def MotifSearchBF ( dna , k , d ) :
    patterns = [ ]
    for i in range ( 0 , len ( dna [ 0 ] ) - k + 1 ) :
        neighbors = bf . Neighbors ( dna [ 0 ][ i : i + k ] , d )
        for j in neighbors:
            count = 0
            for l in dna :
                for i in range ( 0 , len ( l ) - k + 1 ) :
                    if bf . HammingDistance ( j , l [ i : i + k ] ) <= d :
                        count += 1
                        break
            if count == len ( dna ) :
                patterns . append ( j )
    Patterns = [] 
    [ Patterns . append ( x ) for x in patterns if x not in Patterns ] 
    return Patterns

# BRUTE FORCE MOTIF SEARCH: Traverses ALL possible kmers to find the one for 
# which the lowest scoring motif set exists. DOUBLE MINIMISATION PROBLEM
# Input : A list of strings of Dna of equal length, length k of pattern to be found
# INEFFICIENT. 4^k can be narrowed down dramatically. 
# Method may not be improved upon, due to existence of more promising alternative.
def BFMotifSearch ( Dna , k ) : #Median String search
    Motifs = [ ]
    p = [ ]
    d = 10000000
    for i in range ( 4 ** k ) :
        p2 = bf . NumberToPattern( i , k )
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

#Starts off with kmer from first strand. Builds profile and selects best ith motif 
# using the first i-1 motifs. Does so for all possible kmers from first strand
# Returns the best of the motif sets encountered this way
def GreedyMotifSearch ( Dna , k ) :         #O(nt(n+kt)) ~ O ( n^2 t ). << O(n^t)
    best = [ i [ 0 : k ] for i in Dna ]
    bscore = Score ( Count ( best ) )
    score = 0
    t = len ( Dna )
    for i in range ( len ( Dna [ 0 ] ) - k + 1 ) :   # O ( n ) iterations 
        Motif = [ Dna [ 0 ][ i : i + k ] ]
        for i in range ( 1 , t ) :                   # O ( t ) iterations
            profile = Profile ( Motif )              # O ( k t )
            Motif . append ( MostProbableKmer ( Dna [ i ] , profile , k ) [ 0 ] ) # O ( n )
        score = Score ( Count ( Motif ) )
        if ( score < bscore ) :
            bscore = score
            best = Motif
    return best

def GreedyMotifSearchPC ( Dna , k ) :         #O(nt(n+kt)) ~ O ( n^2 t ). << O(n^t)
    best = [ i [ 0 : k ] for i in Dna ]
    bscore = Score ( CountPC ( best ) )
    score = 0
    t = len ( Dna )
    for i in range ( len ( Dna [ 0 ] ) - k + 1 ) :   # O ( n ) iterations 
        Motif = [ Dna [ 0 ][ i : i + k ] ]
        for i in range ( 1 , t ) :                   # O ( t ) iterations
            profile = ProfilePC ( Motif )              # O ( k t )
            Motif . append ( MostProbableKmer ( Dna [ i ] , profile , k ) [ 0 ] ) # O ( n )
        score = Score ( CountPC ( Motif ) )
        if ( score < bscore ) :
            bscore = score
            best = Motif
    return best

def RandomizedMotifSearch ( Dna , k ) : # 1 round doesn't guarantee much
    best = RandomMotifs ( Dna , k )
    bscore = Score ( CountPC ( best ) )
    M = best
    while ( True ) : 
        profile = ProfilePC ( M )
        M = Motifs ( profile , Dna , k )
        score = Score ( CountPC ( M ) )
        if ( score < bscore ) :
            bscore = score
            best = M
        else : 
            return best , bscore
     
def RepeatedRandomizedMotifSearch ( Dna , k , N ) : 
    BestMotifs , bscore = RandomizedMotifSearch ( Dna , k )
    for i in range ( N - 1 ):
        M , score = RandomizedMotifSearch ( Dna , k )
        if bscore > score :
            BestMotifs = M
            bscore = score
        if ( i % 20 == 0 ) : 
            print ( i )
    return BestMotifs , bscore

def GibbsSampler ( Dna , k , N ) :
    best = RandomMotifs ( Dna , k )
    bscore = Score ( CountPC ( best ) )
    score = None
    M = best 
    t = len ( Dna ) 
    for i in range ( N - 1 ) :
        j = random . randint ( 0 , t - 1 )
        M . pop ( j )
        profile = ProfilePC ( M )
        M . insert ( j , ProfileGeneratedString ( Dna [ j ] , profile , k ) ) 
        #M [ j ] = ProfileGeneratedString ( Dna [ j ] , profile , k )
        score = Score ( CountPC ( M ) )
        if ( score < bscore ) :
            best = M
            bscore = score
    return best , bscore

def RepeatedGibbsSampler ( Dna , k , n_inner , N ) :
    best , bscore = GibbsSampler ( Dna , k , n_inner )
    for i in range ( N - 1 ) :
        M , score = GibbsSampler ( Dna , k , n_inner ) 
        if ( bscore > score ) : 
            best = M
            bscore = score
        if ( i % 2 == 0 ) :
            print ( i ) 
    return best , bscore

###############################################################################
#                            TESTING GROUNDS                                  #
###############################################################################
