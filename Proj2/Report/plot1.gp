set   autoscale
set term pdf
set output 'plot1.pdf'
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Approximation Ratio"
set xlabel "n"
set ylabel "ratio"
set yr [1.0:1.5]
set xr [0:210]

stats 'results.dat' using ($3/$2) prefix "A"

plot  "results.dat" using 1:($3/$2) title 'ratio' with points linecolor rgb "blue" , \
 4.0/3.0 title "4/3"  with lines linecolor rgb "red", \
 A_mean title "Mean" with lines linecolor rgb "green"