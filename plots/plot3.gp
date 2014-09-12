set   autoscale
set term pdf
set output 'plot3.pdf'
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Affine Linear Comparison"
set xlabel "n"
set ylabel "Affine/Linear"

plot  "data.data" using 1:($3/$2) title 'Ratio' with linespoints linecolor rgb "blue"