TARGET = matrix
OBJS = main.o matrix.o gauss.o gauss_mod.o relax.o input1.o input2.o hilbert.o

CFLAGS = -std=gnu99 -Iinclude -g
LDFLAGS = -lm

OBJS := $(addprefix build/, $(OBJS))

all: $(TARGET)

$(TARGET): $(OBJS)
	gcc -o $@ $^ $(LDFLAGS)

$(OBJS): | build

build: 
	mkdir build

build/%.o: %.c
	gcc -c -o $@ $^ $(CFLAGS)

clean:
	rm $(TARGET) build -rf
