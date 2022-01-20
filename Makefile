build:
	mpicc src/*.c -o fox -lm -g

time:
	mpicc src/*.c -o fox -lm -g -D TIME