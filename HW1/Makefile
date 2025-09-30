# Makefile for matrix-vector multiplication performance comparison

# Compiler and flags
CXX = g++
CXXFLAGS_NONE = -std=c++11
CXXFLAGS_O2 = -std=c++11 -O2
CXXFLAGS_O3 = -std=c++11 -O3

# Source files
SOURCES = sep_matvec.cpp sep_matvec_unroll.cpp con_matvec.cpp con_matvec_manual_unroll.cpp con_matvec_gcc_unroll.cpp

# Target executables for each optimization level
TARGETS_NONE = sep_matvec_none sep_matvec_unroll_none con_matvec_none con_matvec_manual_unroll_none con_matvec_gcc_unroll_none
TARGETS_O2 = sep_matvec_o2 sep_matvec_unroll_o2 con_matvec_o2 con_matvec_manual_unroll_o2 con_matvec_gcc_unroll_o2
TARGETS_O3 = sep_matvec_o3 sep_matvec_unroll_o3 con_matvec_o3 con_matvec_manual_unroll_o3 con_matvec_gcc_unroll_o3

ALL_TARGETS = $(TARGETS_NONE) $(TARGETS_O2) $(TARGETS_O3)

# Default target
all: $(ALL_TARGETS)

# No optimization builds
sep_matvec_none: sep_matvec.cpp
	$(CXX) $(CXXFLAGS_NONE) $< -o $@

sep_matvec_unroll_none: sep_matvec_unroll.cpp
	$(CXX) $(CXXFLAGS_NONE) $< -o $@

con_matvec_none: con_matvec.cpp
	$(CXX) $(CXXFLAGS_NONE) $< -o $@

con_matvec_manual_unroll_none: con_matvec_manual_unroll.cpp
	$(CXX) $(CXXFLAGS_NONE) $< -o $@

con_matvec_gcc_unroll_none: con_matvec_gcc_unroll.cpp
	$(CXX) $(CXXFLAGS_NONE) $< -o $@

# O2 optimization builds
sep_matvec_o2: sep_matvec.cpp
	$(CXX) $(CXXFLAGS_O2) $< -o $@

sep_matvec_unroll_o2: sep_matvec_unroll.cpp
	$(CXX) $(CXXFLAGS_O2) $< -o $@

con_matvec_o2: con_matvec.cpp
	$(CXX) $(CXXFLAGS_O2) $< -o $@

con_matvec_manual_unroll_o2: con_matvec_manual_unroll.cpp
	$(CXX) $(CXXFLAGS_O2) $< -o $@

con_matvec_gcc_unroll_o2: con_matvec_gcc_unroll.cpp
	$(CXX) $(CXXFLAGS_O2) $< -o $@

# O3 optimization builds
sep_matvec_o3: sep_matvec.cpp
	$(CXX) $(CXXFLAGS_O3) $< -o $@

sep_matvec_unroll_o3: sep_matvec_unroll.cpp
	$(CXX) $(CXXFLAGS_O3) $< -o $@

con_matvec_o3: con_matvec.cpp
	$(CXX) $(CXXFLAGS_O3) $< -o $@

con_matvec_manual_unroll_o3: con_matvec_manual_unroll.cpp
	$(CXX) $(CXXFLAGS_O3) $< -o $@

con_matvec_gcc_unroll_o3: con_matvec_gcc_unroll.cpp
	$(CXX) $(CXXFLAGS_O3) $< -o $@

# Run all benchmarks
run: $(ALL_TARGETS)
	@echo "========================================"
	@echo "Running Matrix-Vector Multiplication Benchmarks"
	@echo "========================================"
	@echo
	@echo "--- NO OPTIMIZATION ---"
	@echo "Separate Storage:"
	@./sep_matvec_none
	@echo
	@echo "Separate Storage with GCC Pragma Unroll:"
	@./sep_matvec_unroll_none
	@echo
	@echo "Contiguous Storage:"
	@./con_matvec_none
	@echo
	@echo "Contiguous Storage with Manual Unroll:"
	@./con_matvec_manual_unroll_none
	@echo
	@echo "Contiguous Storage with GCC Pragma Unroll:"
	@./con_matvec_gcc_unroll_none
	@echo
	@echo "--- O2 OPTIMIZATION ---"
	@echo "Separate Storage:"
	@./sep_matvec_o2
	@echo
	@echo "Separate Storage with GCC Pragma Unroll:"
	@./sep_matvec_unroll_o2
	@echo
	@echo "Contiguous Storage:"
	@./con_matvec_o2
	@echo
	@echo "Contiguous Storage with Manual Unroll:"
	@./con_matvec_manual_unroll_o2
	@echo
	@echo "Contiguous Storage with GCC Pragma Unroll:"
	@./con_matvec_gcc_unroll_o2
	@echo
	@echo "--- O3 OPTIMIZATION ---"
	@echo "Separate Storage:"
	@./sep_matvec_o3
	@echo
	@echo "Separate Storage with GCC Pragma Unroll:"
	@./sep_matvec_unroll_o3
	@echo
	@echo "Contiguous Storage:"
	@./con_matvec_o3
	@echo
	@echo "Contiguous Storage with Manual Unroll:"
	@./con_matvec_manual_unroll_o3
	@echo
	@echo "Contiguous Storage with GCC Pragma Unroll:"
	@./con_matvec_gcc_unroll_o3

# Clean build artifacts
clean:
	rm -f $(ALL_TARGETS)

# Help target
help:
	@echo "Available targets:"
	@echo "  all       - Build all executables with different optimization levels"
	@echo "  run       - Build and run all benchmarks"
	@echo "  clean     - Remove all built executables"
	@echo "  help      - Show this help message"
	@echo
	@echo "Individual targets:"
	@echo "  *_none    - Build with no optimization"
	@echo "  *_o2      - Build with -O2 optimization"
	@echo "  *_o3      - Build with -O3 optimization"

.PHONY: all run clean help