set terminal x11 1 # plots on screen in figure 1
set termoption enhanced
set encoding iso_8859_1 # allows for special symbols like the Angstrom symbol
reset

if (!exists("T")) {
    print "Please enter the temperature on the command line as \'gnuplot -e T=<temperature>\'"
    exit
} else {
    print "T = ".T
}

if (!exists("el")) {
    print "Please enter the element on the command line as \`gnuplot -e \"el=\'<element>\'\"\`"
    exit
} else {
    print "element = ".el
}

if (!exists("r")) {
    print "Please enter the radius of the grain on the command line as \"gnuplot -e r=<radius>\""
    exit
} else {
    print "r = ".r
}

if (!exists("filename")) {
    print "Please enter the filename to plot as \`gnuplot -e \"filename=\'<filename>\'\"\`"
    exit
} else {
    print "filename = ".filename
    axis = system('echo "'.filename.'" | egrep -o "_[0-9][0-9][0-9]_" | cut -c 2-4')
    potential = system('echo "'.filename.'" | egrep -o "_[eb]a[sm][a_](k_)?" | cut -c 2- | head -c -2')
    print "potential = ".potential
}
print "\n"

# Full data plot
full_data_title = el." MSD xy (R = ".r."\305) - T = ".T." K"
fitted_data_title = el." MSD xy (R = ".r."\305) Fitted Region - T = ".T." K"
x_label = 'Time (ps)'
y_label = 'MSD (\305^2)'

set title full_data_title
set xlabel x_label
set ylabel y_label
set yrange [0:]

num_points = system("echo $(($(cat ".filename." | wc -l) - 2))")
plot filename u ($1*0.002):($6) skip 2 notitle

d_x = 0.2

satisfied = 0
set fit logfile "/dev/null"
set fit quiet
while (!satisfied) {
    print "Please click where the linear region begins."
    pause mouse

    if (GPVAL_DATA_X_MIN >= floor(MOUSE_X / d_x * 10) / 10) {
        rmin = 0
    } else {
        rmin = floor(MOUSE_X / d_x)
    }

    print "Now click where the growth stops."
    pause mouse
    #print MOUSE_X, MOUSE_Y

    if (num_points <= ceil(MOUSE_X / d_x * 10) / 10) {
        rmax = num_points
    } else {
        rmax = ceil(MOUSE_X / d_x)
    }

    ftemp(x) = m*x + y0
    print "Start point: ".rmin
    print "End point: ".rmax
    print "Fitting to ".(rmax-rmin)." out of ".num_points." data points."
    fit ftemp(x) filename u ($1*0.002):($6) every ::rmin::rmax via m, y0
    replot [rmin * d_x + GPVAL_DATA_X_MIN:rmax * d_x + GPVAL_DATA_X_MIN] ftemp(x) lt rgb "red" lw 2 notitle
    print "Does this match the average data well? (Y|N) "
    pause mouse keypress

    #print MOUSE_KEY # y is 121, Y is 89
    if (MOUSE_KEY == 121 || MOUSE_KEY == 89) {
        satisfied = 1
    } else {
        satisfied = 0
    }
}
D = m / 4

set print 'diffusion_coeffs.txt' append
print axis, " ", T, " ",r, " ", potential, " ", D
