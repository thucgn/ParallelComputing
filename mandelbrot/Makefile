CC=gcc
CFLAGS=
CLDFLAGS= -lGL -lGLU -lglut

all : main1 main2 main4 main8

main1 : mandelbrot.c
	$(CC) $(CFLAGS) $(CLDFLAGS) -DPRINT_TIME -o $@ $<

main2 : mandelbrot.c
	$(CC) $(CFLAGS) $(CLDFLAGS) -fopenmp  -DPRINT_TIME -DUSE_OPENMP -DTHREADS=2 -o $@ $<

main4 : mandelbrot.c
	$(CC) $(CFLAGS) $(CLDFLAGS) -fopenmp  -DPRINT_TIME -DUSE_OPENMP -DTHREADS=4 -o $@ $<

main8 : mandelbrot.c
	$(CC) $(CFLAGS) $(CLDFLAGS) -fopenmp  -DPRINT_TIME -DUSE_OPENMP -DTHREADS=8 -o $@ $<

clean : 
	rm -f *.o ./main*
