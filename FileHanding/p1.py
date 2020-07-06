##################### For dealing with .txt and other similar file inputs ( custom inputs ) #########################

with open ( '<InsertFileName' , mode = 'r' ) as f :
    f . readline ( 10000000000 ) # Reads one long line

with open ( '/Users/sashank/Desktop/Data/d.txt' , mode = 'r') as f:
    s = f . read ( ) . split ( ) # Read list with one entry per line

with open ( '/Users/sashank/Desktop/Data/out.txt' , mode = 'w' ) as f : #Store list line by line
    for i in l :        
        f . write ( str ( i ) + '\n' )
