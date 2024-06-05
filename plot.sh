#!/bin/bash


gnuplot >/dev/null 2>&1 <<EOF  
set terminal png nocrop enhanced size 800,600 
set output 'PDD.png'
set colorsequence classic
#default|podo|

set xlabel "Profondeur [cm]"
set ylabel "Dose"
set grid

#  10FF
# Prof [mm]	Raw Value
# 100.00		5.8683E+00

#  06FF
# Prof [mm]	Raw Value
# 100.00		3.8492E+00

set style line 1 lt 2 lc 'red'   lw 1 pt 7 ps 2
set style line 2 lt 2 lc 'green' lw 1 pt 7 ps 2
set style line 3 lt 1 lc 'red'   lw 2 pt 3
set style line 4 lt 1 lc 'green' lw 2 pt 3
set style line 5 lt 2 lc 'blue'  lw 1 pt 7 ps 2

Norm10=5.8683E+00
Norm06=3.8492E+00

plot "PDD_06FF.dat" u (\$1/10):(\$2/Norm06) t 'PDD 06FF' w l ls 3, \
     "PDD_10FF.dat" u (\$1/10):(\$2/Norm10) t 'PDD 10FF' w l ls 4

set xrange[8:12]

set output 'PDD_zoomed.png'

# pdd06(x)= a06*x+b06
# pdd10(x)= a10*x+b10
# pdd06(x)= a06*x+1-10*a06
# pdd10(x)= a10*x+1-10*a10

# pdd(x)= A*exp(-k*x)
# A = exp(10k)

pdd06(x)= exp(k06*(10-x))
pdd10(x)= exp(k10*(10-x))

fit pdd06(x) "PDD_06FF.dat" u (\$1/10):(\$2/Norm06) via k06

fit pdd10(x) "PDD_10FF.dat" u (\$1/10):(\$2/Norm10) via k10

# To for to pass through (prof=10cm; 1.0)
b06=1-10*a06
b10=1-10*a10

plot "PDD_06FF.dat" u (\$1/10):(\$2/Norm06) t 'PDD 06FF' w lp ls 1, \
     pdd06(x) t 'Fit 06FF' w l ls 3, \
     "PDD_10FF.dat" u (\$1/10):(\$2/Norm10) t 'PDD 10FF' w lp ls 2, \
     pdd10(x) t 'Fit 10FF' w l ls 4

set output 'PDD_fitted.png'
plot pdd06(x) t 'PDD Fit 06FF' w l ls 3, \
     pdd10(x) t 'PDD Fit 10FF' w l ls 4

EOF

#qiv PDD.png &
#qiv PDD_zoomed.png &
#qiv PDD_fitted.png &

mv fit.log fitPDD.log

gnuplot >/dev/null 2>&1 <<EOF
#set terminal png nocrop enhanced size 800,600 
set terminal pngcairo enhanced size 800,600 
set output 'Profil06FF.png'
set colorsequence classic
#default|podo|
set xlabel "Distance [cm]"
set ylabel "Dose"
set grid

set style line 1 lt 2 lc 'red'   lw 1 pt 7 ps 2
set style line 2 lt 2 lc 'green' lw 1 pt 7 ps 2
set style line 3 lt 1 lc 'red'   lw 2 pt 3
set style line 4 lt 1 lc 'green' lw 2 pt 3
set style line 5 lt 2 lc 'red'   lw 1 dt 2 pt 7 ps 2
set style line 6 lt 2 lc 'green' lw 1 dt 2 pt 7 ps 2
set style line 7 lt 1 lc 'red'   lw 2 dt 2 pt 3
set style line 8 lt 1 lc 'green' lw 2 dt 2 pt 3

#  10FF
# Dist [mm]	Raw Value
#  0.00		5.8683E+00  X
#  0.00		5.8533E+00  Y

#  06FF
# Dist [mm]	Raw Value
#  0.00		3.8537E+00  X 
#  0.00		3.8537E+00  Y

Norm10=(5.8683E+00 + 5.8533E+00)/2

Norm06=3.8537E+00

set output 'Profil06FF.png'
plot "PRX_06FF.dat" u (\$1/10):(\$2/Norm06) t 'ProfilX 06FF' w l ls 3, \
     "PRY_06FF.dat" u (\$1/10):(\$2/Norm06) t 'ProfilY 06FF' w l ls 7

set output 'Profil10FF.png'
plot "PRX_10FF.dat" u (\$1/10):(\$2/Norm10) t 'ProfilX 10FF' w l ls 4, \
     "PRY_10FF.dat" u (\$1/10):(\$2/Norm10) t 'ProfilY 10FF' w l ls 8


set xrange[-2:2]

set output 'Profil06FFZoomed.png'
plot "PRX_06FF.dat" u (\$1/10):(\$2/Norm06) t 'ProfilX 06FF' w l ls 3, \
     "PRY_06FF.dat" u (\$1/10):(\$2/Norm06) t 'ProfilY 06FF' w l ls 7

set output 'Profil10FFZoomed.png'
plot "PRX_10FF.dat" u (\$1/10):(\$2/Norm10) t 'ProfilX 10FF' w l ls 4 , \
     "PRY_10FF.dat" u (\$1/10):(\$2/Norm10) t 'ProfilY 10FF' w l ls 8


pro06x(x)=1-qx06*x*x
pro06y(x)=1-qy06*x*x
pro06(x)=1-q06*x*x

pro10x(x)=1-qx10*x*x
pro10y(x)=1-qy10*x*x
pro10(x)=1-q10*x*x

q06=(qx06+qy06)/2
q10=(qx10+qy10)/2

fit pro06x(x) "PRX_06FF.dat" u (\$1/10):(\$2/Norm06) via qx06
fit pro06y(x) "PRY_06FF.dat" u (\$1/10):(\$2/Norm06) via qy06
q06=(qx06+qy06)/2

fit pro10x(x) "PRX_10FF.dat" u (\$1/10):(\$2/Norm10) via qx10
fit pro10y(x) "PRY_10FF.dat" u (\$1/10):(\$2/Norm10) via qy10
q10=(qx10+qy10)/2

set output 'Profil06FFFitted.png'
plot "PRX_06FF.dat" u (\$1/10):(\$2/Norm06) t 'ProfilX 06FF' w lp ls 1, \
     "PRY_06FF.dat" u (\$1/10):(\$2/Norm06) t 'ProfilY 06FF' w lp ls 5, \
     pro06(x) t 'Fit' w l ls 3

set output 'Profil10FFFitted.png'
plot "PRX_10FF.dat" u (\$1/10):(\$2/Norm10) t 'ProfilX 10FF' w lp ls 2, \
     "PRY_10FF.dat" u (\$1/10):(\$2/Norm10) t 'ProfilY 10FF' w lp ls 6, \
     pro10(x) t 'Fit' w l ls 4
EOF

mv fit.log fitProfil.log

#qiv Profil06FF.png &
#qiv Profil10FF.png &
#qiv Profil06FFZoomed.png &
#qiv Profil10FFZoomed.png &
#qiv Profil06FFFitted.png &
#qiv Profil10FFFitted.png &

echo
echo "Résultats des fits"
echo "=================="

cat fitPDD.log | grep "+/-"
cat fitProfil.log | grep "+/-"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate py310

./determination_k_FFF.py -d -n 80 > results080.log
cat results080.log | grep ^* | tr -d '*' > points_cyl_080.dat

#>/dev/null 2>&1 <<EOF
gnuplot <<EOF
set view equal xyz
splot "points_cyl_080.dat" u 1:2:3 w d
pause 20
EOF

./determination_k_FFF.py -n 160 > results.log

echo
echo "Résultats des facteurs 6FFF et 10FFF"
echo "===================================="

cat results.log | grep FFF
echo
echo
