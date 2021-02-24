set terminal pngcairo
set termoption enhanced
set encoding iso_8859_1
set output "MSD_plot.png"
reset

set title "MSD"
set xlabel "Time (ps)"
set ylabel "MSD (\305^2)"
set yrange [0:35]

if (!exists("file")) {
    print "Please enter the file to plot as \`gnuplot -e \"file=\'<filename>\'\"\`"
    exit
}

if (!exists("dt")) {
    print "Assuming dt=0.002"
    dt = 0.002
}

set key top left
plot 1/0 with points pt 7 ps 0.5 lc rgb "blue" title "x", file u ($1*dt):($2) skip 2 with points pointtype 0 lc rgb "blue" notitle, \
    1/0 with points pt 7 ps 0.5 lc rgb "red" title "y", file u ($1*dt):($3) skip 2 with points pointtype 0 lc rgb "red" notitle, \
    1/0 with points pt 7 ps 0.5 lc rgb "green" title "z", file u ($1*dt):($4) skip 2 with points pointtype 0 lc rgb "green" notitle
