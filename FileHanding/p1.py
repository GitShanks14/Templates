import time

##################### For dealing with .txt ad other similar file inputs ( custom inputs ) #########################

t0 = time . time ( )
with open ( '<InsertFileName' , mode = 'r' ) as f :
    # Readline ( upper_limit ) returns line by line
    # Perform operations on the data read to avoid exceptions. Do it outside if data = Null works
    print ( "Time Elapsed : " , ( time . time ( ) - t0 ) )
    t0 = time . time ( )

#Perform intermediate operations and repeat the same in write mode
