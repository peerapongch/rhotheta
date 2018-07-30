REM: Compile file for Microsoft C 
cl  /Fegenetree sublike.c uniform.c multiple.c filearg.c integrat.c two.c  general.c lu.c paths.c seq2tr.c treepic.c
cl /Fetreepic -DTREEPIC -DGENETREE -DPSVIEW treepic.c
cl /Feseq2tr -DSEQTOTR seq2tr.c
del *.obj
