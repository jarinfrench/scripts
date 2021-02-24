set terminal x11 1 enhanced
set encoding iso_8859_1
reset

if (!exists("file")) {
    print "Please enter the filename to plot as \`gnuplot -e \"file=\'<filename>\'\"\`"
    exit
}

if (!exists("dt")) {
    print "Assuming dt = 0.002"
    dt = 0.002
}

set xlabel 'Time (ps)'
set ylabel 'MSD (\305^2)'
set yrange [0:]

d_x = system("head -n4 ".file." | tail -n2 | awk '{s=$0; getline; print $0-s}'") * dt

x_using = "($1*dt):2"
y_using = "($1*dt):3"
z_using = "($1*dt):4"
xy_using = "($1*dt):($2+$3)"
xz_using = "($1*dt):($2+$4)"
yz_using = "($1*dt):($3+$4)"
xyz_using = "($1*dt):5"

num_points = system("echo $(($(cat ".file." | wc -l) - 2))")
set key top left
plot 1/0 w p pt 7 lc rgb "#FF0000" title "x",\
     1/0 w p pt 7 lc rgb "#00FF00" title "y",\
     1/0 w p pt 7 lc rgb "#0000FF" title "z",\
     1/0 w p pt 7 lc rgb "#FFFF00" title "xy",\
     1/0 w p pt 7 lc rgb "#FF00FF" title "xz",\
     1/0 w p pt 7 lc rgb "#00FFFF" title "yz",\
     1/0 w p pt 7 lc rgb "#FFFFFF" title "xyz",\
     file u ($1*dt):2 w p pt 0 lc rgb "#FF0000" notitle,\
     "" u ($1*dt):3 w p pt 0 lc rgb "#00FF00" notitle,\
     "" u ($1*dt):4 w p pt 0 lc rgb "#0000FF" notitle,\
     "" u ($1*dt):($2+$3) w p pt 0 lc rgb "#FFFF00" notitle,\
     "" u ($1*dt):($2+$4) w p pt 0 lc rgb "#FF00FF" notitle,\
     "" u ($1*dt):($3+$4) w p pt 0 lc rgb "#00FFFF" notitle,\
     "" u ($1*dt):5 w p pt 0 lc rgb "#FFFFFF" notitle


set fit logfile "/dev/null"
set fit quiet

diffusion_calcs = "x_using y_using z_using xy_using xz_using yz_using xyz_using"

iter = 0
do for [diff in diffusion_calcs] {
    satisfied = 0
    set terminal x11 2 enhanced
    unset key
    plot file u value(diff) w p pt 0 lc rgb "#FFFFFF" notitle

    while (!satisfied) {
        print "Click the start point for calculating the diffusion coefficient for".diff
        pause mouse

        if (GPVAL_DATA_X_MIN >= floor(MOUSE_X / d_x * 10) / 10) {
            rmin = 0
        } else {
            rmin = floor(MOUSE_X / d_x)
        }

        print "Click the end point"
        pause mouse

        if (num_points <= ceil(MOUSE_X / d_x * 10) / 10) {
            rmax = num_points
        } else {
            rmax = ceil(MOUSE_X / d_x)
        }

        print "Start point: ".rmin
        print "End point: ".rmax
        fit line_eq(x) file u value(diff) every ::rmin::rmax via a1, a0
        replot [rmin*d_x + GPVAL_DATA_X_MIN:rmax*dx+GPVAL_DATA_X_MIN] line_eq(x) lt rgb "red" lw 2 notitle
        print "Is this fit satisfactory for calculating the diffusion coefficient?"
        pause mouse keypress
        if (MOUSE_KEY == 121 || MOUSE_KEY == 89) {
            satisfied = 1
        } else {
            satisfied = 0
        }
        iter = iter + 1
    }
    print (iter % 3 + 1)
    D = a1 / (iter % 3 + 1)

    print diff." ".D
}
