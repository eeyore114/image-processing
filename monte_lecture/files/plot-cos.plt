set terminal pdfcairo
set output 'plot-cos.pdf'
set xrange [-2*pi:2*pi]
set yrange [-1.5:1.5]
set samples 500
plot cos(x)
