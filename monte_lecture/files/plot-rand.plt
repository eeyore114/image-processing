set terminal pdfcairo
set output 'plot-rand.pdf'
set xrange [-2:2]
set yrange [-2:2]
plot "uniform_rand_hist.txt" using 1:2

