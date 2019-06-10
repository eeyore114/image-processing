set terminal pdfcairo
set output 'plot-sinc.pdf'
set xrange [-10*pi:10*pi]
set yrange [-1:1.2]
set samples 500
plot sin(x)/x
