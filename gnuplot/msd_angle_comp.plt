set terminal pngcairo size 900, 1600 font "Gill Sans,14"
set termoption enhanced
set encoding iso_8859_1

set xrange [0:2500]
ymaxdefault = 10000000

nangles = system("find . -name '20degree' -o -name '30degree' -o -name '45degree' -o -name 'sigma7' | wc -l")
ntemps = system("find . -name 'T*' | cut -d '/' -f 3 | sort | uniq | wc -l")
ndirs = system("find . -name 'dir_*' | cut -d '/' -f 5 | sort | uniq | wc -l")
nfiles = system("find . -name 'MSD_U.dat' | wc -l")

angles = system("find . -name '20degree' -o -name '30degree' -o -name '45degree' -o -name 'sigma7' | cut -d '/' -f 2 | sort")
temps = system("find . -name 'T*' | cut -d '/' -f 3 | sort | uniq")
dirs = system("find . -name 'dir_*' | cut -d '/' -f 5 | sort | uniq")

print "Found ".nfiles." files in ".nangles." angles and ".ntemps." temperatures over ".ndirs." runs"

if (!exists("dt")) {
    print "Assuming dt=0.002"
    dt = 0.002
}

do for [dir in dirs] {
    set output "MSD_angle_comp_".dir.".png"
    set multiplot title "MSD" font "Gill Sans, 36" layout ntemps,nangles margins 0.15, 0.95, 0.15, 0.95 spacing 0.04, 0.02
    row = 1
    col = 1
    do for [T in temps] {
        set yrange[0:ymaxdefault]
        ymaxrange = 0
        do for [angle in angles] {
            file = angle."/".T."/large_r/".dir."/MSD_U.dat"
            if (!file_exists(file)) {
                continue
            }
            stats file u 2 name "X_".angle nooutput
            stats file u 3 name "Y_".angle nooutput
            stats file u 4 name "Z_".angle nooutput
            ymaxrange = max(ymaxrange, max(value("X_".angle."_max"), max(value("Y_".angle."_max"), value("Z_".angle."_max"))))
        }

        yticsetting = getincr(ymaxrange, 4, 1)
        xticsetting = getincr(2500, 5, 2)
        set yrange[0:ymaxrange]

        do for [angle in angles] {
            file = angle."/".T."/large_r/".dir."/MSD_U.dat"

            if (col == 1) {
                if (!file_exists(file)) {
                    set ytics yticsetting format "%g" nomirror
                } else {
                    set ytics yticsetting format "%g"
                }
                set label T at graph -0.4,0.5 center rotate font "Gill Sans, 18"
            } else {
                if (!file_exists(file)) {
                    set ytics yticsetting format "" nomirror
                } else {
                    set ytics yticsetting format ""
                }
                unset label
            }
            if (row == ntemps) {
                if (!file_exists(file)) {
                    set xtics yticsetting format "%g" nomirror rotate by -60 left
                } else {
                    set xtics yticsetting format "%g" rotate by -60 left
                }
            } else {
                if (!file_exists(file)) {
                    set xtics yticsetting format "" nomirror
                } else {
                    set xtics yticsetting format ""
                }
                unset xlabel
            }
            if (col == 1 && row == ceil(ntemps/2.0)) {
                set ylabel "MSD (\305^2)" font "Gill Sans,22"
            } else {
                unset ylabel
            }
            if (row == ntemps && col == ceil(nangles/2.0)) {
                set xlabel "Time (ps)" font "Gill Sans,22" center
            } else {
                unset xlabel
            }
            if (row == 1) {
                set title angle offset 0,-1 font "Gill Sans, 18"
            } else {
                unset title
            }
            if (!file_exists(file)) {
                unset border
                if (row == ntemps && col == 1) {
                    set key left bmargin
                    plot 1/0 w p pt 7 lc rgb "#FF0000" title "x",\
                         1/0 w p pt 7 lc rgb "#00FF00" title "y",\
                         1/0 w p pt 7 lc rgb "#0000FF" title "z"
                } else {
                    plot NaN notitle
                }
                set border
                col = col + 1
                continue
            }
            if (row == ntemps && col == 1) {
                set key left bmargin
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
        }
        row = row + 1
        col = 1
    }
    unset multiplot
}
