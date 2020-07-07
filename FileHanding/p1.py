##################### For dealing with .txt and other similar file inputs ( custom inputs ) #########################

with open ( '<InsertFileName' , mode = 'r' ) as f :
    f . readline ( 10000000000 ) # Reads one long line

with open ( '/Users/sashank/Desktop/Data/d.txt' , mode = 'r') as f:
    s = f . read ( ) . split ( ) # Read list with one entry per line

with open ( '/Users/sashank/Desktop/Data/out.txt' , mode = 'w' ) as f : #Store list line by line
    for i in l :        
        f . write ( str ( i ) + '\n' )

        
###################################################################################################################
# RANDOM ONE TIME USE STUFF FOR A COURSE THAT I DON'T WANT TO DELETE BECAUSE I MIGHT HAVE TO REUSE SOON           #
###################################################################################################################
# WRITE GRAPH IN FILE
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
 
f . write ( '->' . join ( [ str ( val ) for val in Path ] ) ) # WRITE PATH IN FILE
