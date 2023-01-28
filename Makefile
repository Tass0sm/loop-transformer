
CC=gcc
BASEOPTS=-I. -O3 polybench.c -D POLYBENCH_TIME -Wno-unknown-pragmas
SEQOPTS=$(BASEOPTS) -D OPT
PAROPTS=$(SEQOPTS) -fopenmp


base:
	$(CC) $(BASEOPTS) $(TARGET).c -o $(TARGET)-base.exe

seq:
	$(CC) $(SEQOPTS) $(TARGET).c -o $(TARGET)-seq.exe

par:
	$(CC) $(PAROPTS) $(TARGET).c -o $(TARGET)-par.exe

clean:
	rm *.exe
