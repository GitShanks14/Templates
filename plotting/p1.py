
import numpy as np
import matplotlib.pyplot as plt

#################################### 1. SAME X AXIS BUT DIFF Y AXES SCALES ###########################################

# Initialise data
x = np . linspace ( 0 , 2 * np . pi , 10000)
y = np . sin ( x ) 
y2 = np . sinh ( x )

# Set spacing
x_ticks = np . linspace ( 0 , 7 , 11 )
y_ticks1 = np . linspace ( -1 , 1 , 11 )
y_ticks2 = np . linspace ( 0 , 300 , 21 ) 

# Create subplots, edit separately and "superposition" the plots
figure , ax1 = plt . subplots ( )
ax2 = ax1 . twinx ( )
curve1, = ax1 . plot ( x , y , color = 'c' , label = 'Sine' )
curve2, = ax2 . plot ( x , y2, color = 'm' , label = 'Sinh' )
curves = [ curve1 , curve2 ]

# Label the graph
ax2 . legend ( curves , [ curve . get_label ( ) for curve in curves ] )
ax1 . set_xlabel ( "Angle" , color = curve2 . get_color ( ) ) 
ax1 . set_ylabel ( "Sine" , color = curve1 . get_color ( ) )
ax2 . set_ylabel ( "Sinh" , color = curve2 . get_color ( ) )

# Graph sheet settings
ax1 . tick_params ( axis = 'x' , color = 'c' )
ax1 . set_xticks ( x_ticks )

'''
# 1. If the desired grid is wrto y axis 1
ax1 . grid ( color = 'g' )

# 2. If the desired grid is wrto y axis 2
ax1 . grid ( color = 'g' )  # To print grid from x axis too
ax2 . grid ( color = 'g' )  # To print grid wrto y axis of plot 2
ax1 . yaxis . grid ( False )# To hide the grid wrto y axis of plot 1 as goal is 
                            # to print grid wrto only y axis of plot 2
'''

plt . title ( "Real vs Imaginary inputs")








################################################# Displaying / Exporting plot ################################################
plt . show ( )
