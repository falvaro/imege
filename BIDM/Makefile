FLAGS=-lm -Wno-unused-result -O3

bidm: libpbm.o libpgm.o bidm.o pgmtextfea.o
	gcc -o bidm bidm.o pgmtextfea.o libpgm.o libpbm.o $(FLAGS)

bidm.o: bidm.c pgmtextfea.o
	gcc -c bidm.c $(FLAGS)

pgmtextfea.o: pgmtextfea.h pgmtextfea.c libpgm.o
	gcc -c pgmtextfea.c $(FLAGS)

libpgm.o: libpgm.h libpgm.c
	gcc -c libpgm.c $(FLAGS)

libpbm.o: libpbm.h libpbm.c
	gcc -c libpbm.c $(FLAGS)

clean:
	rm -f *.o
