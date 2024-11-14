# Compiler and flags
CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors

# Target executable and object files
TARGET = symnmf
OBJS = symnmf.o

# Default target: build everything
all: $(TARGET) build_ext

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

# Compile symnmf.c into an object file
symnmf.o: symnmf.c symnmf.h
	$(CC) $(CFLAGS) -c symnmf.c

# Build Python extension using setup.py
build_ext:
	python3 setup.py build_ext --inplace

# Clean generated files
clean:
	rm -f $(TARGET) $(OBJS) *.so
	rm -rf build

# Phony targets to avoid conflicts with files of the same name
.PHONY: all clean build_ext
