##################### For dealing with .txt and other similar file inputs ( custom inputs ) #########################

with open ( '<InsertFileName' , mode = 'r' ) as f :
    f . readline ( 10000000000 ) # Reads one long line

with open ( '/Users/sashank/Desktop/Data/out.txt' , mode = 'w' ) as f :
    for i in l :                       #Writes list entries line by line
        f . write ( str ( i ) + '\n' )

