# Makefile

BINDIR := bin
LIBDIR := lib

CCFLAGS := -pedantic

#CC := g++ -std=c++17
CC:= mpicxx

# src/ (declaration of functions, classes + code)
# main/ (main)
# bin/ (temporary, .o, .exe)
# lib/ (libraries)
# Figures (plots and figures)
# Data (saved data)

# making library
# - static: .a
# - shareable: .so

VPATH = main:src

SRC := $(wildcard src/*.C)
OBJ := $(patsubst %.C, $(BINDIR)/%.o, $(notdir $(SRC)))
INC := $(wildcard src/*.h)
PY := $(wildcard pyt/*.py)
File:= $(wildcard Data/*.txt) $(wildcard Data/*.dat)
Fig:= $(wildcard Figures/*.png) $(wildcard Figures/*.pdf) $(wildcard */*.png)
lib: $(LIBDIR)/libFC.a

$(LIBDIR)/libFC.a: $(OBJ) 
	@echo make lib...
	ar ruv $@ $^
	ranlib $@

%.exe: $(BINDIR)/%.o $(LIBDIR)/libFC.a
	@echo compilink and linking... 
	$(CC) -I src $< -o $(BINDIR)/$@ -L lib -l FC

$(BINDIR)/%.o: %.C | $(INC)
	@echo compiling... $<
	$(CC) -I src -c $< -o $@

python:
	python3 $(PY)

######### clean

tilde := $(wildcard */*~) $(wildcard *~)
exe := $(wildcard */*.exe) $(wildcard *.exe)
obj := $(wildcard */*.o) $(wildcard *.o) $(wildcard */*.so) $(wildcard */*.pcm) $(wildcard */*.d) 

clean:
	@echo cleaning dir...
	rm -f $(exe) $(obj) $(tilde)

cleanFIG:
	@echo cleaning dir...
	rm -f $(File) $(Fig) 
