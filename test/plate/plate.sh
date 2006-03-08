shapes square | transform --sy 0.06251 | transform --tx .031249 | transform --ty -.015 > square.gts
if gerris2D $1 | awk '{ if ($9 < 10.) exit (1); }'; then :
    exit 1
else
    exit 0
fi
