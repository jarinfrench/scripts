save ARG1."_commands.plt"
set terminal push
set terminal png ARG2
set output ARG1

replot
set output
set terminal pop
