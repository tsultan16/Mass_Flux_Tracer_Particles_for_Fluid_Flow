set terminal png


TOP=0.90
DY = 0.23
tSteps=100
tSkip=2
nx=100


set key font "2"

set ytics font "Verdana,5" 

set autoscale

do for [i=0:tSteps-1] { 
  j=i*tSkip
  print "Time Steps Completed =".i
  filename2="Frames/1d_t=".j.".png"
  set output filename2
  set title "Time Step = #".(j)
  set xlabel "x"
  set ylabel "# Tracers"
  plot "output1d.txt" every ::j*nx+1::nx+(j*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle
  unset xlabel
  unset ylabel
  unset title

  unset output
}



