set terminal x11 1
set termoption enhanced
set encoding iso_8859_1
reset

x_label = 'Time (ps)'
y_label = 'Area (\305^2)'
set xlabel x_label
set ylabel y_label
set yrange [0:]

if (axis == 100 || axis == 110) {
    plot_index = 2
} else {
    if (axis == 111) {
        plot_index = 3
    }
}

plot file u 1:plot_index

print "Does the plot index need to be switched? (Y|N) "
call "confirm.plt"

if (confirm) {
    if (plot_index == 2) {
        plot_index = 3
    } else {
        plot_index = 2
    }
    plot file u 1:plot_index
}

cmd1 = sprintf("awk '$1~/^[^#]/{print $%d; exit}' %s",plot_index,file) # gets the first non-commented line, and extracts the <plot_index> column
cmd2 = sprintf("awk 'END{print $%d}' %s",plot_index, file) # gets the last line, and extracts the <plot_index> column
y0 = system(cmd1)
yf = system(cmd2)

cmd3 = sprintf("echo 'scale=4; (%s-%s)/%s*100' | bc",yf, y0, y0) # calculates the percent growth
result = system(cmd3)
if (abs(result) < 25)  {
    cmd4 = sprintf('echo "Approximate shrinkage percentage: \033[01;31m%.2f%\033[0m"',abs(result))
    system(cmd4)
    # print "Approximate shrinkage percentage: ".result."%"
    print "WARNING: not enough growth!"
} else {
    cmd4 = sprintf('echo "Approximate shrinkage percentage: \033[01;32m%.2f%\033[0m"',abs(result))
    system(cmd4)
}

d_x = ceil(0.002*10000)
num_points = system("echo $(($(cat ".file." | wc -l) - 1))")

confirm = 0
set fit logfile "/dev/null"
set fit errorvariables
set fit quiet
while (!confirm) {
    print "Please click where the growth starts becoming linear."
    pause mouse

    offset = GPVAL_DATA_X_MIN / d_x
    if (GPVAL_DATA_X_MIN / d_x - offset >= floor(MOUSE_X / d_x) - offset) {
        rmin = 0
    } else {
        rmin = floor(MOUSE_X / d_x)
    }

    print "Now click where the growth stops."
    pause mouse

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
    fit ftemp(x) file u 1:plot_index every ::rmin::rmax via m, y0

    replot [rmin * 0.002 * 10000:rmax * 0.002 * 10000] ftemp(x) lt rgb "red" lw 2 notitle
    print "Does this match the data well? (Y|N) "

    call "confirm.plt"
}

stats file u 1:plot_index every ::rmin::rmax nooutput
set print "fitted_values.txt" append
out_line = sprintf("%s dA/dt = %8.6f error = %8.6f r_sq = %6.4f fit to points ".rmin.":".rmax, file, m, m_err, STATS_correlation*STATS_correlation)
print out_line

#call "export.plt" "tmp.png"
