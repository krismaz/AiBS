set   autoscale
set term pdf
set output 'plot1.pdf'
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Time is O(n^2)"
set xlabel "n"
set ylabel "time(ms)/n^2"

plot  "data.data" using 1:($2/($1*$1)) title 'Linear' with linespoints linecolor rgb "blue", \
		"data.data" using 1:($3/($1*$1)) title 'Affine' with linespoints linecolor rgb "red"