
extra=

opts=-O3 -Wall -shared -fPIC

saiga12c.so: saiga12c.o SFMT.o

	gcc -Wall -O2 -shared -fPIC saiga12c.o SFMT.o -o saiga12c.so

saiga12c.o: saiga12c.c ccode/kobandersen.c ccode/fredricksonandersen.c\
	ccode/birolimezard.c
	gcc ${extra} ${opts} -c saiga12c.c

SFMT.o: SFMT.c SFMT.h
	gcc ${opts} -c -DMEXP=19937 -include SFMT-params.h SFMT.c


