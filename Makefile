
# Ouput options

#OPT   +=  -DSTELLARAGE
#OPT   +=  -DMETALS
#OPT   +=  -DOUTPUTPOTENTIAL
#OPT   +=  -OUTPUTACCELERATION
#OPT   +=  -DOUTPUTCHANGEOFENTROPY
#OPT   +=  -DOUTPUTTIMESTEP

# simulation options

#OPT   +=  -DLONGIDS
#OPT   +=  -DCOOLING	
#OPT   +=  -DSFR   #only for runs with P-Gadget2

# Additional routines

OPT   +=  -DREAD_GADGET1
OPT   +=  -DTRANSLATIONS_ROTATIOS
#OPT   +=  -DBINS
OPT   +=  -DMEDLEY

CC=gcc
UBICATIONGSL=/home/lquiroga/local
#CFLAGS=-g -I. -I$(UBICATIONGSL)/include $(OPT)
CFLAGS= -I. -O3 -Wall -I$(UBICATIONGSL)/include $(OPT) 
LFLAGS= -lm -L$(UBICATIONGSL)/lib -lgsl -lgslcblas

clean:
	rm -rf *.o* *~ *.out

%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@
	cp $@ /home/lquiroga/local/bin



