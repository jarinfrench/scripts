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
    potential = system('echo "'.basename.'" | egrep -o "_[eb]a[sm][a_](k_)?" | cut -c 2- | head -c -2')
    if (el eq "UO2") {
        print "potential = ".potential
    }
}
print "\n"

# Full data plot
force_vel_data_title = el." Velocity vs Force (r_0 = ".r."\305) - T = ".T." K"
x_label = 'Force (GPa)'
y_label = 'Velocity (m/s)'

set title force_vel_data_title
set xlabel x_label
set ylabel y_label
set key inside bottom right

average_file = system("ls ".basename."_average.txt")

num_list = system("ls -v ".basename."_[0-9]*.txt | egrep -o [0-9][0-9]?\.txt$ | cut -d. -f1")
plot for [i in num_list] basename."_".i.".txt" index 0 u 4:5 title "Run ".i, \
basename."_average.txt" index 0 u 4:5 title "Average"

call "export.plt" basename."_force_vs_velocity.png"
replot
