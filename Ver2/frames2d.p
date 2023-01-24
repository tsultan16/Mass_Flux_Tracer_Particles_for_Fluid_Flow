set terminal png
set pm3d map
set palette color

set dgrid3d 100,100


TOP=0.90
DY = 0.23
tSteps=50
tSkip=1

do for [i=1:tSteps] { 
  j=i*tSkip
  print "Time Steps Completed =".i
  filename2="Frames/frame_t=".j.".png"
  set output filename2

  filename="t=".j.".txt"
  set title "Time Step = #".j
  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:3 
  unset title
  unset xlabel
  unset ylabel
  
  unset output
}



