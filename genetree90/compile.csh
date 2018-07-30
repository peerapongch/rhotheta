# Version 1.0 13/8/97
cc  -O -o genetree sublike.c uniform.c multiple.c filearg.c integrat.c two.c  general.c lu.c paths.c seq2tr.c treepic.c -lm
rm *.o
cc -O -o treepic -DTREEPIC -DGENETREE -DPSVIEW treepic.c -lm
cc -O -o seq2tr -DSEQTOTR seq2tr.c -lm
rm *.o
