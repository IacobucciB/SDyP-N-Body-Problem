CC = gcc
CC = gcc
CFLAGS = -O1 -Wall -std=c99 -Wno-missing-braces -Iinc

LDFLAGS = -Llib -lraylib -lopengl32 -lgdi32 -lwinmm -lpthread

STATIC_LDFLAGS = -Llib -lraylib -lopengl32 -lgdi32 -lwinmm -lpthread \
                 -static -static-libgcc -static-libstdc++

SRC = src/main.c src/nbody.c
OUT = build/main.exe
STATIC_OUT = build/main_static.exe

all:
	$(CC) -o $(OUT) $(SRC) $(CFLAGS) $(LDFLAGS)

static:
	$(CC) -o $(STATIC_OUT) $(SRC) $(CFLAGS) $(STATIC_LDFLAGS)

clean:
	del /Q build\*.exe 2> NUL || true

.PHONY: all clean static