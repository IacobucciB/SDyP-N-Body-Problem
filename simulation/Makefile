CC = gcc
CFLAGS = -O1 -Wall -std=c99 -Wno-missing-braces -Iinc
LDFLAGS = -Llib -lraylib -lopengl32 -lgdi32 -lwinmm -lpthread

SRC = src/main.c src/nbody.c
OUT = build/main.exe

all:
	$(CC) -o $(OUT) $(SRC) $(CFLAGS) $(LDFLAGS)

clean:
	del /Q build\*.exe 2> NUL || true

.PHONY: all clean