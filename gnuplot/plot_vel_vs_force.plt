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

# if (!exists("basename")) {
#     print "Please enter the base file name of the files to plot as \`gnuplot -e \"basename=\'<base file name>\'\"\`"
#     exit
# } else {
#     print "basename = ".basename
#     axis = system('echo "'.basename.'" | egrep -o "_[0-9][0-9][0-9]_" | cut -c 2-4')
#     print "axis = ".axis
#     potential = system('echo "'.basename.'" | egrep -o "_[eb]a[sm][a_](k_)?" | cut -c 2- | head -c -2')
#     if (el eq "UO2") {
#         print "potential = ".potential
#     }
# }
print "\n"

# Full data plot
force_vel_data_title = el." Velocity vs Force (r_0 = ".r."\305) - T = ".T." K"
x_label = 'Inverse radius (\305^{-1})'
y_label = 'Velocity (m/s)'

set title force_vel_data_title
set xlabel x_label
set ylabel y_label
set cblabel "Time (ps)"
set key inside top right

stats "force_velocity_data.txt" using 1 nooutput
set xrange [0:0.1]
set yrange [-1:15]

save "force_vs_velocity_commands.plt"
#set terminal push
set terminal png
set output "force_vs_velocity.png"

plot "force_velocity_data.txt" index 0 u (1/$3):($1 < STATS_max - 120 ? -1*$5 :1/0):1 every :::0::0 with points palette notitle
#replot
#set output
#set terminal pop
#call "export.plt" "force_vs_velocity.png"
#replot
