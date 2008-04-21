
#extra="-DFFid=$i"

#for neighlist in "-Dneighlist" "" ; do
neighlist=""

opts="-O3 -Wall -shared -fPIC $neighlist"

gcc $* $extra $opts -c saiga12c.c
gcc $opts -c -DMEXP=19937 -include SFMT-params.h SFMT.c

# With icc instead of gcc:
#icc -Wall -O2 -shared -fPIC -c dragunov_c.c
#icc -Wall -O2 -shared -fPIC -c -DMEXP=19937 -include SFMT-params.h SFMT.c

# link it together (gcc needed regardless of what compiler you use first)
#gcc -Wall -O2 -shared -fPIC saiga12c.o SFMT.o -o saiga12c.so
gcc -Wall -O2 -shared -fPIC saiga12c.o SFMT.o -o saiga12c$neighlist.so

#done

