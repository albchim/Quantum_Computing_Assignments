
set terminal png size 800,600
set fit quiet


N=1999*100
h_file = 'out_h.txt'
d_file = 'out_d.txt'


P(x) = a*(x**alpha) * exp(-b*(x**beta))

set ylabel 'P(s)'
set xlabel 'Normalized spacings'

# HIST
binwidth = 0.2
set boxwidth binwidth
set style fill solid 0.8
bin(x, width) = width*floor(x/width) + width/2


set table "hist_h.txt"
plot h_file u (bin($1,binwidth)):(1.0/(binwidth*N)) smooth freq with boxes
unset table

# FIT
a = 2
b = 2
alpha = 2
beta = 2 
fit P(x) "hist_h.txt" u 1:2 via a,alpha,b,beta

# PLOT
set output 'herm.png'
set yrange [0:1.02]
set xrange [-0.5:8]
set title "Hermitian matrix"
plot h_file u (bin($1, binwidth)):(1.0/(binwidth*N)) smooth freq w boxes lc rgb 'orange' title 'Nomalized spacings', P(x) with lines lw 3 lc rgb 'royalblue' title 'P(s)'



set table "hist_d.txt"
plot d_file u (bin($1,binwidth)):(1.0/(binwidth*N)) smooth freq with boxes
unset table

# FIT
a = 2
b = 2
alpha = 2
beta = 2 
fit P(x) "hist_d.txt" u 1:2 via a,alpha,b,beta

# PLOT
set output 'diag.png'
set yrange [0:1.02]
set xrange [-0.5:8]
set title "Real Diagonal matrix"
plot d_file u (bin($1, binwidth)):(1.0/(binwidth*N)) smooth freq w boxes lc rgb 'orange' title 'Nomalized spacings', P(x) with lines lw 3 lc rgb 'royalblue' title 'P(s)'

