# set terminal pngcairo font "Gill Sans,9" linewidth 4 rounded fontscale 1.0
set terminal x11
set encoding iso_8859_1
set termoption enhanced

# Line styles: pick pleasing colors, rather than strictly primary colors or
# hard-to-see colors like gnuplot's default yellow. Make the lines thick so they
# are easy to see in small plots on papers
set style line 1 lt rgb "#A00000" lw 2 pt 1
set style line 2 lt rgb "#00A000" lw 2 pt 6
set style line 3 lt rgb "#5060D0" lw 2 pt 2
set style line 4 lt rgb "#F25900" lw 2 pt 9

dx_file = "dx_diffusion_data.txt"
dy_file = "dy_diffusion_data.txt"
dz_file = "dz_diffusion_data.txt"
dxy_file = "dxy_diffusion_data.txt"
dxz_file = "dxz_diffusion_data.txt"
dyz_file = "dyz_diffusion_data.txt"
dxyz_file = "dxyz_diffusion_data.txt"
file_list = "dx_diffusion_data.txt dy_diffusion_data.txt dz_diffusion_data.txt dxy_diffusion_data.txt dxz_diffusion_data.txt dyz_diffusion_data.txt dxyz_diffusion_data.txt"
title_list = "D_x D_y D_z D_{xy} D_{xz} D_{yz} D_{xyz}"
section_list = "Dx Dy Dz Dxy Dxz Dyz Dxyz"

sequence_list = "100\\\_20degree 100\\\_30degree 100\\\_45degree 110\\\_20degree 110\\\_30degree 110\\\_45degree 111\\\_20degree 111\\\_30degree 111\\\_45degree 111\\\_Sigma7"
name_list = "100_20degree 100_30degree 100_45degree 110_20degree 110_30degree 110_45degree 111_20degree 111_30degree 111_45degree 111_Sigma7"
axes = "100 110 111"
angs = "20degree 30degree 45degree"
seq_100 = "100\\\_20degree 100\\\_30degree 100\\\_45degree" # start index at 0
seq_110 = "110\\\_20degree 110\\\_30degree 110\\\_45degree" # start index at 3
seq_111 = "111\\\_20degree 111\\\_30degree 111\\\_45degree 111\\\_sigma7" # start index at 6
div_list = "ADP\\\_100 ADP\\\_110 ADP\\\_111 MEAM\\\_100 MEAM\\\_110 MEAM\\\_111 Ternary\\\_EAM\\\_100 Ternary\\\_EAM\\\_110 Ternary\\\_EAM\\\_111"

xs = 2 # index for the xs
ys = 5 # index for the ys
yerrs = 6 # index for the y error

# ADP plots
set yrange [-9.5:-1]
set print "diffusion_parameters.txt"
print "potential_axis_angle D_0(m^2/s) q(eV)"
print "ADP"
do for [i=1:words(file_list)] {
    print word(section_list,i)
    do for [j=1:words(name_list)] {
        var_name = "ADP_".word(name_list,j)."_".word(section_list,i)
        stats "adp_".word(file_list,i) every :::j-1::j-1 u xs:ys name var_name nooutput
        intercept = var_name."_intercept"
        slope = var_name."_slope"
        d0 = gprintf("%.2se%.2S", exp(value(intercept))*1e-8) # in units of m^2/s
        q = gprintf("%.2f", value(slope)*-1)
        if (var_name eq "ADP_111_Sigma7") {
            if (strlen(d0) < 8) {
                print var_name."\t\t".d0."\t\t".q
            } else {
                print var_name."\t\t".d0."\t".q
            }
        } else {
            if (strlen(d0) < 8) {
                print var_name."\t".d0."\t\t".q
            } else {
                print var_name."\t".d0."\t".q
            }
        }
    }
    print " "
}

plot for [i=1:words(sequence_list)] "adp_".dx_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i)
call "export.plt" "adp_dx_diffusion.png"
plot for [i=0:words(sequence_list)] "adp_".dy_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "adp_dy_diffusion.png"
plot for [i=0:words(sequence_list)] "adp_".dz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "adp_dz_diffusion.png"
plot for [i=0:words(sequence_list)] "adp_".dxy_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "adp_dxy_diffusion.png"
plot for [i=0:words(sequence_list)] "adp_".dxz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "adp_dxz_diffusion.png"
plot for [i=0:words(sequence_list)] "adp_".dyz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "adp_dyz_diffusion.png"
plot for [i=0:words(sequence_list)] "adp_".dxyz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "adp_dxyz_diffusion.png"

plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::0::0 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_100_20degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::1::1 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_100_30degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::2::2 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_100_45degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::3::3 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_110_20degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::4::4 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_110_30degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::5::5 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_110_45degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::6::6 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_111_20degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::7::7 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_111_30degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::8::8 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_111_45degree_diffusion.png"
plot for [i=1:words(file_list)] "adp_".word(file_list,i) every :::9::9 u xs:ys:yerrs w yerrorbars t "ADP ".word(title_list,i)
call "export.plt" "adp_d_111_sigma7_diffusion.png"

# MEAM plots
print "MEAM"
do for [i=1:words(file_list)] {
    print word(section_list,i)
    do for [j=1:words(name_list)] {
        var_name = "MEAM_".word(name_list,j)."_".word(section_list,i)
        stats "meam_".word(file_list,i) every :::j-1::j-1 u xs:ys name var_name nooutput
        intercept = var_name."_intercept"
        slope = var_name."_slope"
        d0 = gprintf("%.2se%.2S", exp(value(intercept)))
        q = gprintf("%.2f", value(slope)*-1)
        if (var_name eq "MEAM_111_Sigma7") {
            if (strlen(d0) < 8) {
                print var_name."\t\t".d0."\t\t".q
            } else {
                print var_name."\t\t".d0."\t".q
            }
        } else {
            if (strlen(d0) < 8) {
                print var_name."\t".d0."\t\t".q
            } else {
                print var_name."\t".d0."\t".q
            }
        }
    }
    print ""
}

plot for [i=0:words(sequence_list)] "meam_".dx_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "meam_dx_diff.png"
plot for [i=0:words(sequence_list)] "meam_".dy_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "meam_dy_diff.png"
plot for [i=0:words(sequence_list)] "meam_".dz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "meam_dz_diff.png"
plot for [i=0:words(sequence_list)] "meam_".dxy_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "meam_dxy_diff.png"
plot for [i=0:words(sequence_list)] "meam_".dxz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "meam_dxz_diff.png"
plot for [i=0:words(sequence_list)] "meam_".dyz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "meam_dyz_diff.png"
plot for [i=0:words(sequence_list)] "meam_".dxyz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "meam_dxyz_diff.png"

plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::0::0 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_100_20degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::1::1 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_100_30degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::2::2 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_100_45degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::3::3 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_110_20degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::4::4 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_110_30degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::5::5 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_110_45degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::6::6 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_111_20degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::7::7 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_111_30degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::8::8 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_111_45degree_diffusion.png"
plot for [i=1:words(file_list)] "meam_".word(file_list,i) every :::9::9 u xs:ys:yerrs w yerrorbars t "MEAM ".word(title_list,i)
call "export.plt" "meam_d_111_sigma7_diffusion.png"

# Ternary EAM plots
print "Ternary EAM"
do for [i=1:words(file_list)] {
    print word(section_list,i)
    do for [j=1:words(name_list)] {
        var_name = "Ternary_EAM_".word(name_list,j)."_".word(section_list,i)
        stats "ternary_eam_".word(file_list,i) every :::j-1::j-1 u xs:ys name var_name nooutput
        intercept = var_name."_intercept"
        slope = var_name."_slope"
        d0 = gprintf("%.2se%.2S", exp(value(intercept)))
        q = gprintf("%.2f", value(slope)*-1)
        if (var_name eq "Ternary_EAM_111_Sigma7") {
            if (strlen(d0) < 8) {
                print var_name."\t\t".d0."\t\t".q
            } else {
                print var_name."\t\t".d0."\t".q
            }
        } else {
            if (strlen(d0) < 8) {
                print var_name."\t".d0."\t\t".q
            } else {
                print var_name."\t".d0."\t".q
            }
        }
    }
    print ""
}

plot for [i=0:words(sequence_list)] "ternary_eam_".dx_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "ternary_eam_dx_diff.png"
plot for [i=0:words(sequence_list)] "ternary_eam_".dy_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "ternary_eam_dy_diff.png"
plot for [i=0:words(sequence_list)] "ternary_eam_".dz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "ternary_eam_dz_diff.png"
plot for [i=0:words(sequence_list)] "ternary_eam_".dxy_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "ternary_eam_dxy_diff.png"
plot for [i=0:words(sequence_list)] "ternary_eam_".dxz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "ternary_eam_dxz_diff.png"
plot for [i=0:words(sequence_list)] "ternary_eam_".dyz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "ternary_eam_dyz_diff.png"
plot for [i=0:words(sequence_list)] "ternary_eam_".dxyz_file every :::i-1::i-1 u xs:ys:yerrs w yerrorbars t word(sequence_list,i+1)
call "export.plt" "ternary_eam_dxyz_diff.png"

plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::0::0 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_100_20degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::1::1 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_100_30degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::2::2 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_100_45degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::3::3 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_110_20degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::4::4 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_110_30degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::5::5 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_110_45degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::6::6 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_111_20degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::7::7 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_111_30degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::8::8 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_111_45degree_diffusion.png"
plot for [i=1:words(file_list)] "ternary_eam_".word(file_list,i) every :::9::9 u xs:ys:yerrs w yerrorbars t "Ternary EAM ".word(title_list,i)
call "export.plt" "ternary_eam_d_111_sigma7_diffusion.png"

set print "diffusion_parameters_by_potential.txt"
print "Potential_axis_angle\tDx_int\tDx_q\tDy_int\tDy_q\tDz_int\tDz_q\tDxy_int\tDxy_q\tDxz_int\tDxz_q\tDyz_int\tDyz_q\tDxyz_int\tDxyz_q"
potentials = "ADP MEAM Ternary_EAM"
do for [a=1:words(potentials)] {
    do for [i=1:words(name_list)] {
        line = word(potentials, a)."_".word(name_list,i)
        do for [j=1:words(section_list)] {
            var_name = word(potentials,a)."_".word(name_list,i)."_".word(section_list,j)
            intercept = var_name."_intercept"
            slope = var_name."_slope"
            d0 = gprintf("%.2se%.2S", exp(value(intercept)))
            q = gprintf("%.2f", value(slope)*-1)
            line = line."\t".d0."\t".q."\t"
        }
        print line
    }
    print ""
}

titles = titles = system("awk 'NR > 1 !/^#/ {gsub(\"_\",\"-\",$1); print $1}' diffusion_parameters_by_potential.txt")
potentials = "ADP MEAM Ternary\\\_EAM"
set logscale x
# Shows the spread of each diffusion coefficient across all potentials
plot for [i=1:words(div_list)] 'diffusion_parameters_by_potential.txt' u 2:3 every :::i-1::i-1 skip 1 pt 6 lc rgb '#E41A1C' notitle,\
     for [i=1:words(div_list)] 'diffusion_parameters_by_potential.txt' u 4:5 every :::i-1::i-1 skip 1 pt 7 lc rgb '#377EB8' notitle,\
     for [i=1:words(div_list)] 'diffusion_parameters_by_potential.txt' u 6:7 every :::i-1::i-1 skip 1 pt 8 lc rgb '#4DAF4A' notitle,\
     for [i=1:words(div_list)] 'diffusion_parameters_by_potential.txt' u 8:9 every :::i-1::i-1 skip 1 pt 9 lc rgb '#984EA3' notitle,\
     for [i=1:words(div_list)] 'diffusion_parameters_by_potential.txt' u 10:11 every :::i-1::i-1 skip 1 pt 10 lc rgb '#FF7F00' notitle,\
     for [i=1:words(div_list)] 'diffusion_parameters_by_potential.txt' u 12:13 every :::i-1::i-1 skip 1 pt 11 lc rgb '#FFFF33' notitle,\
     for [i=1:words(div_list)] 'diffusion_parameters_by_potential.txt' u 14:15 every :::i-1::i-1 skip 1 pt 13 lc rgb '#A65628' notitle,\
     1/0 with points pt 6 lc rgb '#E41A1C' title "Dx",\
     1/0 with points pt 7 lc rgb '#377EB8' title "Dy",\
     1/0 with points pt 8 lc rgb '#4DAF4A' title "Dz",\
     1/0 with points pt 9 lc rgb '#984EA3' title "Dxy",\
     1/0 with points pt 10 lc rgb '#FF7F00' title "Dxz",\
     1/0 with points pt 11 lc rgb '#FFFF33' title "Dyz",\
     1/0 with points pt 13 lc rgb '#A65628' title "Dxyz"
call "export.plt" "diffusion_coeff_spread.png"
