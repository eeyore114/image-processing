set terminal pdfcairo
set output 'rand.pdf'
set xrange [-1:1]
set yrange [-1:1]
plot sin(2x)

