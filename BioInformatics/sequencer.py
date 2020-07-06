
################################################################################
#                       GENOME SEQUENCER                                       #
################################################################################

import numpy as np
import random 
import sys

sys . path . insert ( 0 , '/Users/sashank/Desktop/Python/Libraries/Bioinformatics/BioUtils' )

import BioUtils as bf


# Input  : A genome
# Output : Genome -> list of all kmers contained in genome 
def Composition ( Pattern , k ) :
    c = [ ] 
    for i in range ( len ( Pattern ) - k + 1 ) :
        c . append ( Pattern [ i : i + k ] )
    return c

# Input  : A Genome
# Output : List of all kmers contained in genome w no repition 
def UniqueComposition ( Pattern , k ) :
    c = [ ] 
    for i in range ( len ( Pattern ) - k + 1 ) :
        if Pattern [ i : i + k ] not in c :
            c . append ( Pattern [ i : i + k ] )
    return c

# Input  : A Path in a graph, represented as a list of nodes of k-mers
# Output : The genome corresponding to the path
def PathToGenome ( Path ) :
    Genome = Path [ 0 ] 
    for i in range ( 1 , len ( Path ) ) :
        Genome += Path [ i ][ -1 :  ]
    return Genome

# Input  : Kmer 
# Output : First k-1 bases
def Prefix ( Genome ) :
    return Genome [ : -1 ]

# Input  : Kmer 
# Output : Last k-1 bases
def Suffix ( Genome ) :
    return Genome [ 1 : ]
 
# Input  : Reads of a genome
# Output : DiGraph that connects every read to every read it can be connected to
def Graph ( Reads ) :
    AdjList = { }
    for i in Reads :                                # n steps
        if i not in AdjList : 
            AdjList [ i ] = [ ] 
        for j in Reads :                            # n steps
            if ( Suffix ( i ) == Prefix ( j ) ) :   # == -> k steps
                AdjList [ i ] . append ( j ) 
    return AdjList

# Input  : Reads of a genome, used as the edges of the DeBrujin Graph
# Output : DiGraph w nodes = (k-1)-mers and the path -> genome
# O(kn^2)
def DeBrujin ( Edges ) :        # Edges = composition of the Genome. 
    Nodes = [ ]
    for i in Edges :            # O(n) nodes -> O(n) edges here
        p = Prefix ( i )
        s = Suffix ( i )
        if p not in Nodes :     # O(n) search twice, O(k) comparison
            Nodes . append ( p )
        if s not in Nodes : 
            Nodes . append ( s )
    n = len ( Edges )
    AdjList = { }
    for i in range ( len ( Nodes ) ) :  # O(n) steps
        AdjList [ Nodes [ i ] ] = [ ] 
        for j in range ( n ) :          # O(n) steps
            if Prefix ( Edges [ j ] ) == Nodes [ i ] : # O(k) comparison
                AdjList [ Nodes [ i ] ] . append ( Suffix ( Edges [ j ] ) ) 
    return AdjList

################################################################################
#                       TESTING GROUNDS                                        #
################################################################################

import time
t0 = time . time ( )

def timer ( ):
    global t0
    print ( "Execution of block took {} s" . format ( time . time ( ) - t0 ) )
    t0 = time . time ( )



print ( "Execution started" )


with open ( '/Users/sashank/Desktop/Data/d.txt' , mode = 'r') as f:
    Edges = f . read ( ) . split ( ) # Read list with one entry per line


timer ( )    
AdjList = DeBrujin ( Edges )
timer ( )

l = list ( AdjList . items ( ) )
l . sort ( )

with open ( '/Users/sashank/Desktop/Data/out.txt' , mode = 'w' ) as f : #Store list line by line
    for key , val in l :        
        k = len ( val )
        if ( k == 0 ) :
            continue 
        i = 0
        f . write ( str ( key ) + ' -> ' )
        while ( i < k ) :
            f . write ( val [ i ] )
            if ( i != k - 1 ) :
                f . write ( ',' )
            i += 1
        f . write ( '\n' )
timer ( )
    
    
    
    
    
    