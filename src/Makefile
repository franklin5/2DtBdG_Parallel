CC=mpicxx
eigenbin = /home/ld7/bin	
CFLAGS=-c -Wall -I${eigenbin}
LDFLAGS=-limf -lm
SRC=main.cpp grst.cpp gapnumberintegral.cpp tBdG.cpp lgwt.cpp
OBJ=$(SRC:.cpp=.o)
EXE=tbdg

all: $(SRC) $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)
.cpp.o: 
	$(CC) $(CFLAGS) $< -o $@

doc:
	doxygen Doxyfile

touch: 
	touch *.cpp *.h
clean:
	rm *.o $(EXE) *~ *#

