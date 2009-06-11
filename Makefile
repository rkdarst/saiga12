
extra=

opts=-O3 -Wall -shared -fPIC

saiga12c.so: saiga12c.o SFMT.o
	gcc -Wall -O2 -shared -fPIC ${extra} saiga12c.o SFMT.o -o saiga12c.so

saiga12c.o: saiga12c.c ccode/kobandersen.c ccode/fredricksonandersen.c\
	ccode/birolimezard.c ccode/east.c ccode/energy_bm.c \
	ccode/ctcc.c ccode/ctccclassic.c ccode/energy_ctcc.c
	gcc ${extra} ${opts} -c saiga12c.c

SFMT.o: SFMT.c SFMT.h
	gcc ${opts} -c -DMEXP=19937 -include SFMT-params.h SFMT.c

test:
	python tests/unittests_run.py

# add this line to _darcs/prefs/prefs to have darcs auto-test on record:
# test make darcs-test
darcs-test: saiga12c.so
	mkdir run-dir 
	ln -s .. run-dir/saiga12 
	PYTHONPATH=${PWD}/run-dir/ python tests/unittests_run.py

