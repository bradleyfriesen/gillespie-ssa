MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Define color codes
GREEN = '\033[0;32m'
RED = '\033[0;31m'
YELLOW = '\033[0;33m'
RESET = '\033[0m'

# Remove
RM := rm -f

# Directories
SRC := src
OBJ := obj

# Search all source directories
VPATH := $(shell find $(SRC) -type d)

# Compiler and flags
FC := gfortran
FFLAGS := -O3 -fopenmp
LD := $(FC)
LDFLAGS := -O3 -fopenmp
LDLIBS :=

# Sources
MODS := \
	kind_parameters.f90 \
	core/command_line.f90 \
	core/model_interface.f90 \
	core/model_factory.f90 \
	core/SSA.f90 \
	core/statistics.f90 \
	models/birth_death/birth_death_model.f90 \
	models/homodimer/homodimer_model.f90

PROGS := \
	example.f90

MODOBJS := $(addprefix $(OBJ)/,$(notdir $(MODS:.f90=.o)))
PROGOBJS := $(addprefix $(OBJ)/,$(notdir $(PROGS:.f90=.o)))
TARGETS := $(PROGS:.f90=)

.PHONY: all clean
.DEFAULT_GOAL := all

# All targets
all: $(TARGETS)

$(OBJ):
	mkdir -p $(OBJ)

# Rule for compiling sources
$(OBJ)/%.o: %.f90 | $(OBJ)
	@echo -e $(YELLOW)Compiling $@$(RESET)
	$(FC) $(FFLAGS) -J$(OBJ) -c $< -o $@

# Programs depend on all modules
$(PROGOBJS):$(MODOBJS)

# Linking programs
%: $(MODOBJS) $(OBJ)/%.o
	@echo -e $(YELLOW)Linking $@$(RESET)
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Module dependencies
$(OBJ)/kind_parameters.o:
$(OBJ)/command_line.o: \
	$(OBJ)/model_factory.o
$(OBJ)/model_interface.o: \
	$(OBJ)/kind_parameters.o
$(OBJ)/model_factory.o: \
	$(OBJ)/model_interface.o \
	$(OBJ)/birth_death_model.o \
	$(OBJ)/homodimer_model.o
$(OBJ)/statistics.o: \
	$(OBJ)/kind_parameters.o
$(OBJ)/SSA.o: \
	$(OBJ)/kind_parameters.o \
	$(OBJ)/model_interface.o
$(OBJ)/birth_death_model.o: \
	$(OBJ)/kind_parameters.o \
	$(OBJ)/model_interface.o
$(OBJ)/homodimer_model.o: \
	$(OBJ)/kind_parameters.o \
	$(OBJ)/model_interface.o

# Clean up build files
clean:
	$(RM) $(MODOBJS) $(PROGOBJS) $(TARGETS) $(OBJ)/*.mod
