OBJ= modelprob.o obsData.o hmmresults.o hmm.o baumwelch.o hillclimbing.o main.o

all: hmm

hmm: ${OBJ}
	g++ -o hmm ${OBJ} -lpopt -lnlopt
clean:
	rm ${OBJ} hmm

