set   autoscale
set term pdf
set output 'plot3.pdf'
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Time is O(n^2)"
set xlabel "n"
set ylabel "time(ms)/n^2"
set yr [0:0.03]

plot  "data.data" using 1:($3/$2) title 'Linear' with linespoints linecolor rgb "blue"