
CFLAGS=-O3

opts=-Wall -shared -fPIC

saiga12c.so: saiga12c.o SFMT.o
	gcc ${opts} ${CFLAGS} saiga12c.o SFMT.o -o saiga12c.so

saiga12c.o: saiga12c.c ccode/kobandersen.c ccode/fredricksonandersen.c\
	ccode/birolimezard.c ccode/east.c ccode/energy_bm.c \
	ccode/ctcc.c ccode/ctccclassic.c ccode/energy_ctcc.c \
	ccode/spiral.c ccode/energy_squareplaquette.c \
	ccode/spinmontecarlo.c
	gcc ${opts} ${CFLAGS} -c saiga12c.c

SFMT.o: SFMT.c SFMT.h
	gcc ${opts} ${CFLAGS} -DMEXP=19937 -include SFMT-params.h -c SFMT.c

test:
	python tests/unittests_run.py

# add this line to _darcs/prefs/prefs to have darcs auto-test on record:
# test make darcs-test
darcs-test: saiga12c.so
	mkdir run-dir 
	ln -s .. run-dir/saiga12 
	PYTHONPATH=${PWD}/run-dir/ python tests/unittests_run.py




#libcutil.py: _cutil.c _cutil.h _cutil.so
#	python -m ctypeslib.h2xml -c $(PWD)/_cutil.h -o __tmp-mpi.xml
#	python -m ctypeslib.xml2py __tmp-mpi.xml -l $(PWD)/_cutil.so \
#	    -o libcutil.py
#	rm __tmp-mpi.xml
#_cutil.so: _cutil.c
#	gcc ${opts} ${extra} -c _cutil.c
#	gcc ${opts} ${extra} _cutil.o -o _cutil.so
