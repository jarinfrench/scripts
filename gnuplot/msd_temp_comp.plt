set terminal pngcairo size 1920, 1600 font "Gill Sans,14"
set termoption enhanced
set encoding iso_8859_1
set output "MSD_T_comp_plot.png"
reset

set xrange [0:2500]

if (!exists("ymax")) {
    set yrange [0:35]
} else {
    set yrange [0:ymax]
}

ntemps = system("find . -name 'T*' | wc -l")
ndirs = system("find . -name 'dir_*' | cut -d '/' -f 4 | sort | uniq | wc -l")
nfiles = system("find . -name 'MSD_U.dat' | wc -l")

temps = system("find . -name 'T*' | cut -d '/' -f 2 | sort")
dirs = system("find . -name 'dir_*' | cut -d '/' -f 4 | sort | uniq")

print "Found ".nfiles." files in ".ntemps." temperatures over ".ndirs." runs"

if (!exists("dt")) {
    print "Assuming dt=0.002"
    dt = 0.002
}

yticson = 1
xticson = ndirs
set multiplot title "MSD" font "Gill Sans, 36" layout ndirs,ntemps margins 0.15, 0.95, 0.15, 0.95 spacing 0.02, 0.02
unset key
row = 1
col = 1
do for [dir in dirs] {
    do for [T in temps] {
        file = T."/large_r/".dir."/MSD_U.dat"
        if (!file_exists(file)) {
            set multiplot next
            continue
        }
        if (col == 1) {
            set ytics format "%g"
        } else {
            set ytics format ""
        }
        if (row == ndirs) {
            set xtics format "%g" rotate by -45 left

        } else {
            set xtics format ""
            unset xlabel
        }
        if (col == 1 && row == ceil(ndirs/2.0)) {
            set ylabel "MSD (\305^2)" font "Gill Sans,22"
        } else {
            unset ylabel
        }
        if (row == ndirs && col == ceil(ntemps/2.0)) {
            set xlabel "Time (ps)" font "Gill Sans,22" left
        } else {
            unset xlabel
        }
        if (row == 1) {
            set title T offset 0,-1
        } else {
            unset title
        }
        if (row == 1 && col == 1) {
            set key top left
            plot 1/0 w p pt 7 lc rgb "#FF0000" title "x", file u ($1*dt):2 w p pt 0 lc rgb "#FF0000" notitle,\
                1/0 w p pt 7 lc rgb "#00FF00" title "y", file u ($1*dt):3 w p pt 0 lc rgb "#00FF00" notitle,\
                1/0 w p pt 7 lc rgb "#0000FF" title "z", file u ($1*dt):4 w p pt 0 lc rgb "#0000FF" notitle
        } else {
            unset key
            plot file u ($1*dt):2 w p pt 0 lc rgb "#FF0000" notitle,\
                "" u ($1*dt):3 w p pt 0 lc rgb "#00FF00" notitle,\
                "" u ($1*dt):4 w p pt 0 lc rgb "#0000FF" notitle
        }

        col = col + 1
        if (col > ntemps) {
            row = row + 1
            col = 1
        }
    }
}
unset multiplot

set terminal ""
set terminal pngcairo font "Gill Sans,8" size 1800,400
set termoption enhanced
set encoding iso_8859_1
set output "MSD_T_comp_plot_dir1.png"

set multiplot title "MSD" font "Gill Sans,10" layout 2,ntemps margins 0.05, 0.975, 0.15, 0.95 spacing 0.02, 0.02
n = 1
do for [i=1:ntemps] {
    set multiplot next
}
do for [T in temps] {
    file = T."/large_r/dir_1/MSD_U.dat"
    if (!file_exists(file)) {
        print "Could not find MSD file in dir_1 of ".T."/large_r"
        break
    }
    set xtics format "%g" rotate by -45 left
    set title T offset 0,-1
    set tmargin at screen 0.9
    if (n == 1) {
        set key top left
        set ytics format "%g"
        set ylabel "MSD (\305^2)"
        plot 1/0 w p pt 7 ps 0.5 lc rgb "#FF0000" title "x", file u ($1*dt):2 w p pt 0 lc rgb "#FF0000" notitle,\
            1/0 w p pt 7 ps 0.5 lc rgb "#00FF00" title "y", file u ($1*dt):3 w p pt 0 lc rgb "#00FF00" notitle,\
            1/0 w p pt 7 ps 0.5 lc rgb "#0000FF" title "z", file u ($1*dt):4 w p pt 0 lc rgb "#0000FF" notitle
    } else {
        unset key
        unset ylabel
        set ytics format ""
        if (n == ceil(ntemps/2.0)) {
            set xlabel "Time (ps)"
        } else {
            unset xlabel
        }
        plot file u ($1*dt):2 w p pt 0 lc rgb "#FF0000" notitle,\
            "" u ($1*dt):3 w p pt 0 lc rgb "#00FF00" notitle,\
            "" u ($1*dt):4 w p pt 0 lc rgb "#0000FF" notitle
    }
    n = n + 1
}
unset multiplot
