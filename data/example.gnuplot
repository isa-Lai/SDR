# example.gnuplut : configuration for plotting (change as needed)

reset                                   # reset
set size ratio 1                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 2,2 scale 1.0,1.0  # set two plots for this figure

# re-sampled mono
set ylabel 'Spectrum (dB/Hz)'            # set y-axis label
set xlabel 'Frequency (Hz)'             # set x-axis label
set yrange [-1.5:1.5]                       # set y plot range
set xrange [0:1000]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/fm_demod.dat' using 1:2 with lines lt 1 lw 1 lc rgb '#880000' notitle

# re-sampled mono
set ylabel 'Spectrum (dB/Hz)'            # set y-axis label
set xlabel 'Frequency (Hz)'             # set x-axis label
set yrange [-1:1]                       # set y plot range
set xrange [0:1000]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/fm_demod_q.dat' using 1:2 with lines lt 1 lw 1 lc rgb '#880000' notitle


# freq domain (PSD)
set ylabel 'Spectrum (dB/Hz)'            # set y-axis label
set xlabel 'Frequency (KHz)'             # set x-axis label
set yrange [-2:2]                       # set y plot range
set xrange [-2:2]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/demod_time.dat' using 1:2 with points lt 1 lw 1 lc rgb '#880000' notitle


unset multiplot