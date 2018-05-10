# label_loop
# This loop adds a label to a plot by pressing the left mouse key.
# If you are not convinced with your chosen position, just click the mouse key
# again and it will be positioned at another place. If you are finished, just
# press another key.
#
# AUTHOR: Hagen Wierstorf
# UPDATED: 19 July 2017 by Jarin French ($0 deprecated)

# Initialize a label number
if (!exists("label_number")) label_number = 1;

# Waiting for the key press
pause mouse any ARG1

# Check if the left mouse key is pressed and add the given label to the plot.
# Otherwise stop the loop and count the added label
if( MOUSE_BUTTON==1 ) \
    set label label_number ARG1 at MOUSE_X,MOUSE_Y center;\
    print "\n ".ARG1." at ",MOUSE_X,MOUSE_Y;\
    replot;\
    reread;\
else \
    label_number = label_number+1;\
    print "\n"
