
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

def EulerianCycle ( Graph , NewStart = 0 ) :
    MainCycle = [ ]
    Cycle = [ NewStart ]
    node = NewStart 
    while Graph [ node ] :                  # Walk
        node = Graph [ node ] . pop ( 0 )
        Cycle . append ( node )
    for node in Cycle :
        if ( Graph [ node ] ) :
            MainCycle . extend ( EulerianCycle ( Graph , node ) )
        else :
            MainCycle . append ( node )
    return MainCycle

def EulerianPath ( Graph ) :
    count = { }
    elen = 0
    for key , val in Graph . items ( ) :
        if key not in count :
            count [ key ] = -len ( val )
        else :
            count [ key ] -= len ( val )
        elen += len ( val )
        for i in val : 
            if i not in count :
                count [ i ] = 1
            else : 
                count [ i ] += 1
    start = None
    end = None
    for key , val in count . items ( ) :
        if val != 0 :
            if val == -1 and start == None :
                start = key
            elif val == 1 and end == None : 
                end = key
            else : 
                print ( "ERROR: NO EULERIAN PATH PRESENT" )
                return
    if end not in Graph :
        Graph [ end ] = [ ]
    Graph [ end ] . append ( start )
    Path = EulerianCycle ( Graph , start )
    pos = None
    for i in range ( len ( Path ) - 1 ) :
        if Path [ i ] == end and Path [ i + 1 ] == start:
            pos = i + 1
            break
    i = 1 
    n = len ( Path )
    FPath = [ ]
    while i < n :
        FPath . append ( Path [ ( pos + i ) % n ] )
        i += 1
    return FPath

def StringReconstruction ( Patterns , k ):
    db = DeBrujin( Patterns )
    DubGraph = { }
    for key , val in db . items ( ) :
        k2 = PatternToNumber ( key ) 
        DubGraph [ k2 ] = [ ]
        for i in val :
            DubGraph [ k2 ] . append ( PatternToNumber ( i ) )
    RawPath = EulerianPath ( DubGraph )
    Path = [ ]
    for i in RawPath : 
        Path . append ( NumberToPattern ( i , k - 1 ) )
    Genome = PathToGenome ( Path )
    return Genome
        
def KUniversalCircularString ( k ) :
    Edges = [ ]
    for i in range ( 2 ** k ) :
        b = bin ( i )
        b = b [ 2 : ]
        while len ( b ) < k :
            b = '0' + b
        Edges . append ( b )
    db = DeBrujin ( Edges ) 
    DubGraph = { }
    for key , val in db . items ( ) :
        k2 = int ( key , 2 )
        DubGraph [ k2 ] = [ ]
        for i in val :
            DubGraph [ k2 ] . append ( int ( i , 2 ) )
    RawPath = EulerianCycle ( DubGraph , 0 )
    Path = [ ]
    for i in RawPath : 
        b = bin ( i ) 
        b = b [ 2 : ]
        while ( len ( b ) < k - 1 ) :
            b = '0' + b 
        Path . append ( b )
    string = PathToGenome ( Path )
    string = string [ 0 : len ( string ) - k + 1 ]
    return string

################################################################################
#                       TESTING GROUNDS                                        #
################################################################################
