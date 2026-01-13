# Default compiler and flags
CXX = g++
CXXFLAGS = -O2 -std=c++17

# Set these via environment or directly here
PYTHIA8_CONFIG := pythia8-config
HEPMC3_INC := -I$(HEPMC3_DIR)/include
HEPMC3_LIB := -L$(HEPMC3_DIR)/lib -lHepMC3

# Get all .cc files in the directory
SOURCES := $(wildcard *.cc)
TARGETS := $(basename $(SOURCES))

# Default rule: build all executables
all: $(TARGETS)

# Pattern rule to build each executable
%: %.cc
	$(CXX) $(CXXFLAGS) `$(PYTHIA8_CONFIG) --cxxflags --libs` \
	$(HEPMC3_INC) $(HEPMC3_LIB) \
	-Wl,-rpath,$(HEPMC3_DIR)/lib \
	-o $@ $<

# Clean
clean:
	rm -f $(TARGETS)

