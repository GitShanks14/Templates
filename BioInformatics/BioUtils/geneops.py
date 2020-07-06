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