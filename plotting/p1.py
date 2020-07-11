
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
ax1 . set_yticks ( y_ticks1 )
ax2 . set_yticks ( y_ticks2 ) 

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

################################################ Different X axis scale but same Y axis scale ###############################

# Initialise data
x1 = np . sin ( y )
x2 = np . sinh ( y ) 
y = np . linspace ( 0 , 7 , 15 )

# Set spacing
y_ticks = np . linspace ( 0 , 7 , 15 )
x_ticks1 = np . linspace ( -1 , 1 , 11 )
x_ticks2 = np . linspace ( 0 , 300 , 7 ) 

# Create subplots, edit separately and "superposition" the plots
figure , ax1 = plt . subplots ( )
ax2 = ax1 . twiny ( )
curve1, = ax1 . plot ( x1 , y , color = 'c' , label = 'Sine Inverse' )
curve2, = ax2 . plot ( x2 , y , color = 'm' , label = 'Sinh Inverse' )
curves = [ curve1 , curve2 ]

# Label the graph
ax2 . legend ( curves , [ curve . get_label ( )  for curve in curves ] )
ax1 . set_ylabel ( "Angle" , color = curve2 . get_color ( ) ) 
ax1 . set_xlabel ( "ArcSine" , color = curve1 . get_color ( ) )
ax2 . set_xlabel ( "ArcSinh" , color = curve2 . get_color ( ) )

# Graph sheet settings
ax1 . tick_params ( axis = 'x' , color = curve1 . get_color ( ) )
ax2 . tick_params ( axis = 'x' , color = curve2 . get_color ( ) )
ax1 . set_xticks ( x_ticks1 )
ax1 . set_xticks ( x_ticks2 )
ax1 . set_yticks ( y_ticks  )
'''
# 1. If the desired grid is wrto y axis 1
ax1 . grid ( color = 'g' )

# 2. If the desired grid is wrto y axis 2
ax1 . grid ( color = 'g' )  # To print grid from x axis too
ax2 . grid ( color = 'g' )  # To print grid wrto y axis of plot 2
ax1 . xaxis . grid ( False )# To hide the grid wrto y axis of plot 1 as goal is 
                            # to print grid wrto only y axis of plot 2

plt . title ( "Real vs Imaginary inputs")







################################################# Displaying / Exporting plot ################################################
plt . show ( )
