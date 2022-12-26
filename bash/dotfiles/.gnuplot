# Location of helpful gnuplot scripts : Location of helpful configuration files
set loadpath "/home/jarinf/projects/scripts/bin:/home/jarinf/.config/gnuplot"

e = exp(1)
sqrt2 = sqrt(2)
min(a,b) = (a < b) ? a : b
max(a,b) = (a > b) ? a : b
gauss(x,m,s) = exp(-0.5*( (x-m)/s )**2 )/sqrt(2*pi*s**2)

# Helper function to determine if a file exists
file_exists(file) = system("[ -f '".file."' ] && echo '1' || echo '0'") + 0

# Helper functions to handle creating a specific number of tics on a plot
# Taken from https://stackoverflow.com/a/46958680
# It looks like it only really works well for ranges that are factors of 10 or
# 5 * factors of 10 (e.g. 50, 1e-12, 500, 5e20, but not 505, 115, etc))
endsinone(n) = strstrt(gprintf("%g",n), "1")
getincr(range, maxincr, guess) = range/guess < maxincr ? guess : (endsinone(guess) ? getincr(range,maxincr, 5*guess) : getincr(range, maxincr, 2*guess))
line_eq(x) = a1*x + a0
parabola_eq(x) = a2 * x * x + a1 * x + a0
