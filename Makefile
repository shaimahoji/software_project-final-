# Compiler and flags
CC = gcc                                                # The C compiler to use
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors   # Compiler flags for strict C standard and warnings

# Target executable and object files
TARGET = symnmf  # Name of the final executable
OBJS = symnmf.o

# Default target: build everything (executable and Python extension)
all: $(TARGET) build_ext

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm   # With math library

# Compile symnmf.c into an object file
symnmf.o: symnmf.c symnmf.h
	$(CC) $(CFLAGS) -c symnmf.c

# Build Python extension using setup.py
build_ext:
	python3 setup.py build_ext --inplace      # Build the Python C extension in the current directory

# Clean generated files
clean:
	rm -f $(TARGET) $(OBJS) *.so   # Remove the executable, object files, and shared objects
	rm -rf build                   # Remove the build directory created by setup.py

# Phony targets to avoid conflicts with files of the same name
.PHONY: all clean build_ext
