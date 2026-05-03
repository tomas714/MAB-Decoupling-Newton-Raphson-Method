CC     = gcc
CFLAGS = -std=c11 -Wall -Werror -Ofast -g
LDLIBS = -lm
OUTPUT = app

all: $(OUTPUT)

$(OUTPUT): main.o matrix.o mabdec.o
	@echo "Building Executable ..."
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@

%.o: %.c
	@echo "Building object file $@ ..."
	$(CC) $(CFLAGS) -c $< -o $@
	@echo ""

clean:
	@echo "Cleanup ..."
	rm -rf *.o $(OUTPUT)
	@echo ""

rebuild: clean all