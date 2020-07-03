import time

# Define a global variable t0 and store time . time ( ) in it at the start of the program. 
# On calling timer ( ), the execution time between the previous call and the current call will be printed. The first call will print execution time
# from the definition of t0. 
def timer ( ) :
    global t0
    print ( "Execution of block took {} s" . format ( time . time ( ) - t0 ) )
    t0 = time . time ( )
    
   
