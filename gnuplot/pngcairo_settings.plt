set terminal pngcairo font "Gill Sans,9" linewidth 4 rounded fontscale 1.0
set encoding iso_8859_1
set termoption enhanced
# line style for axes
set style line 80 lt rgb "#808080"

# line style for grid
set style line 81 lt 0 # dashed
set style line 81 lt rgb "#808080" # grey

set grid back linestyle 81
set border 3 back linestyle 80 # remove border on top and right, and put it in grey
set xtics nomirror
set ytics nomirror

# Line styles: pick pleasing colors, rather than strictly primary colors or
# hard-to-see colors like gnuplot's default yellow. Make the lines thick so they
# are easy to see in small plots on papers
set style line 1 lt rgb "#A00000" lw 2 pt 1
set style line 2 lt rgb "#00A000" lw 2 pt 6
set style line 3 lt rgb "#5060D0" lw 2 pt 2
set style line 4 lt rgb "#F25900" lw 2 pt 9
