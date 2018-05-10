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

if (!exists("basename")) {
    print "Please enter the base file name of the files to plot as \`gnuplot -e \"basename=\'<base file name>\'\"\`"
    exit
} else {
    print "basename = ".basename
    axis = system('echo "'.basename.'" | egrep -o "_[0-9][0-9][0-9]_" | cut -c 2-4')
    print "axis = ".axis
    if (axis == 100 || axis == 111) {
        plot_index = 3
    } else {
        if (axis == 110) {
            plot_index = 2
        } else {
            print "Unknown axis choice.  Assuming plot_index = 3"
            plot_index = 3
        }
    }
}
print "\n"

# Full data plot
full_data_title = el." Grain Growth (R = ".r."\305) - T = ".T." K"
fitted_data_title = el." Grain Growth (R = ".r."\305) Fitted Region - T = ".T." K"
x_label = 'Time (ps)'
y_label = 'Area (\305^2)'

num_files = system("ls -v ".basename."_[0-9]*.txt | wc -l")
first_file = system("ls -v ".basename."_[0-9]*.txt | head -n 1")
set title full_data_title
set xlabel x_label
set ylabel y_label
set yrange [0:]
if (num_files == 1) {
    print "Only plotting one data set.  Plot will not be saved."
    set terminal x11 1 persist
    plot first_file u 1:plot_index title "Run 1"
    exit
}

average_file = system("ls ".basename."_average.txt")

num_list = system("ls -v ".basename."_[0-9]*.txt | egrep -o [0-9][0-9]?\.txt$ | cut -d. -f1")
num_points = system("echo $(($(cat ".average_file." | wc -l) - 1))")
# At this point, I am just going to assume that the time step is .002 ps, with output every 10000 time steps
# I may try to generalize it later
d_x = ceil(0.002*10000)

plot for [i in num_list] basename."_".i.".txt" u 1:plot_index title "Run ".i, \
basename."_average.txt" u 1:plot_index title "Average"

call "export.plt" basename."_full_data.png"
replot

satisfied = 0
set fit logfile "/dev/null"
set fit quiet
while (!satisfied) {
    print "Please click where the growth starts becoming linear."
    pause mouse

    if (GPVAL_DATA_X_MIN >= floor(MOUSE_X / d_x)) {
        rmin = int(GPVAL_DATA_X_MIN)
    } else {
        rmin = floor(MOUSE_X / d_x)
    }

    print "Now click where the growth stops."
    pause mouse
    #print MOUSE_X, MOUSE_Y

    if (num_points <= ceil(MOUSE_X / d_x)) {
        rmax = num_points
    } else {
        rmax = ceil(MOUSE_X / d_x)
    }

    if (rmin > rmax) {
        continue
    }

    if (rmin < 0) {
        continue
    }
    
    ftemp(x) = m*x + y0
    print "Start point: ".rmin
    print "End point: ".rmax
    print "Fitting to ".(rmax-rmin)." out of ".num_points." data points."
    fit ftemp(x) basename."_average.txt" u 1:plot_index every ::rmin::rmax via m, y0
    replot [rmin * 0.002 * 10000:rmax * 0.002 * 10000] ftemp(x) lt rgb "red" lw 2 notitle
    print "Does this match the average data well? (Y|N) "
    pause mouse keypress

    #print MOUSE_KEY # y is 121, Y is 89
    if (MOUSE_KEY == 121 || MOUSE_KEY == 89) {
        satisfied = 1
    } else {
        satisfied = 0
    }
}


# Now for the fitting
set terminal x11 2

log_name = basename."_fit.log"
set fit logfile log_name
set fit errorvariables
f(x) = dadt * x + b
fit f(x) basename."_average.txt" u 1:plot_index every ::rmin::rmax via dadt, b
set key font ",10"
set xrange [rmin * 0.002 * 10000:rmax * 0.002 * 10000]

plot for [i in num_list] basename."_".i.".txt" u 1:plot_index every ::rmin::rmax title "Run ".i, \
basename."_average.txt" u 1:plot_index every ::rmin::rmax title "Average" lc rgb '#90EE90', \
f(x) lt rgb 'red' lw 2 title 'Fit'

satisfied = 0
while (!satisfied) {
    print "Pick a point to show the slope"
    pause mouse
    set arrow 1 from MOUSE_X,MOUSE_Y to MOUSE_X,f(MOUSE_X) nohead front lc rgb 'black'
    set arrow 2 from MOUSE_X,MOUSE_Y to (MOUSE_Y-b)/dadt,MOUSE_Y nohead front lc rgb 'black'
    replot

    print "Pick where the center of the slope label should go"
    pause mouse
    set label 1 at MOUSE_X,MOUSE_Y sprintf('dA/dt = %3.5g',dadt)." \261 ".sprintf('%3.5g',dadt_err) center front font ",10"
    replot

    print "Is this a good spot? (Y|N)"
    pause mouse keypress

    #print MOUSE_KEY # y is 121, Y is 89
    if (MOUSE_KEY == 121 || MOUSE_KEY == 89) {
        satisfied = 1
    } else {
        satisfied = 0
        #plot for [i=1:num_files] basename."_".i.".txt" u 1:plot_index every ::rmin::rmax title "Run ".i, \
        #basename."_average.txt" u 1:plot_index every ::rmin::rmax title "Average" lc rgb '#90EE90', \
        #f(x) lt rgb 'red' lw 2 title 'Fit'
    }
}

call "export.plt" basename."_fitted_data.png"
