# Makefile for runica_c - Darwin (macOS) platform
# Uses Apple Accelerate framework for BLAS/LAPACK

# Directories
SRCDIR    = src
OBJDIR    = obj

# Compiler and flags
CC        = gcc
CFLAGS    = -O3 -I$(SRCDIR)
LIBS      = -framework Accelerate -lm

# Source files (all .c files in src/)
SRCS      = $(SRCDIR)/runica_simple.c $(SRCDIR)/rtGetInf.c $(SRCDIR)/rtGetNaN.c \
            $(SRCDIR)/rt_nonfinite.c main.c

# Object files (in obj/ directory)
OBJS      = $(OBJDIR)/runica_simple.o $(OBJDIR)/rtGetInf.o $(OBJDIR)/rtGetNaN.o \
            $(OBJDIR)/rt_nonfinite.o $(OBJDIR)/main.o

# Target executable (at root)
TARGET    = runica_darwin

# Default target
all: $(OBJDIR) $(TARGET)

# Create obj directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Link executable
$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(LIBS)

# Compile source files
$(OBJDIR)/runica_simple.o: $(SRCDIR)/runica_simple.c
	$(CC) -c $(SRCDIR)/runica_simple.c -o $(OBJDIR)/runica_simple.o $(CFLAGS)

$(OBJDIR)/rtGetInf.o: $(SRCDIR)/rtGetInf.c
	$(CC) -c $(SRCDIR)/rtGetInf.c -o $(OBJDIR)/rtGetInf.o $(CFLAGS)

$(OBJDIR)/rtGetNaN.o: $(SRCDIR)/rtGetNaN.c
	$(CC) -c $(SRCDIR)/rtGetNaN.c -o $(OBJDIR)/rtGetNaN.o $(CFLAGS)

$(OBJDIR)/rt_nonfinite.o: $(SRCDIR)/rt_nonfinite.c
	$(CC) -c $(SRCDIR)/rt_nonfinite.c -o $(OBJDIR)/rt_nonfinite.o $(CFLAGS)

$(OBJDIR)/main.o: main.c
	$(CC) -c main.c -o $(OBJDIR)/main.o $(CFLAGS)

# Clean build artifacts
clean:
	rm -f $(OBJDIR)/*.o $(TARGET)

# Clean everything including directories
distclean: clean
	rm -rf $(OBJDIR)

.PHONY: all clean distclean
