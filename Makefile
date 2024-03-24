LIB = -L/Users/quanzhou/local/lib  ## please specify the local of GSL
INC = -I/Users/quanzhou/local/include  ## location of GSL
SHELL = /bin/sh
CC = g++
FLAGS = -g -O3
CFLAGS = -lm -lgsl -lgslcblas

TARGET = fastBVSR
SRC = $(shell echo src/*.cpp)
HEADER = $(shell echo src/*.h)
OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(OBJ): %.o: %.cpp
	$(CC) $(FLAGS) -c -o $@ $< $(INC)  

$(TARGET): ${OBJ}
	$(CC) $(FLAGS) -o $@ $^ $(CFLAGS) $(INC) $(LIB)

.PHONY: clean

clean:
	-rm -f fastBVSR $(OBJ)

