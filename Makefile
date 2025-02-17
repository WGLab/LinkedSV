.PHONY: all pigz fermikit clean

# Directories
BIN_DIR   := bin
SRC_DIR   := src
HTSLIB_DIR:= $(SRC_DIR)/htslib-1.3
LIB_DIR   := $(SRC_DIR)/lib
INC_DIR   := $(SRC_DIR)/include

# Compiler and flags
CXX      := g++
CXXFLAGS := -g -O0 -std=c++11 -I $(INC_DIR)

# Default target: build everything
all: pigz fermikit $(BIN_DIR)/cluster_reads $(BIN_DIR)/extract_barcode_info \
     $(BIN_DIR)/output_bam_coreinfo $(BIN_DIR)/remove_sparse_nodes \
     $(BIN_DIR)/cal_hap_read_depth_from_bcd21 $(BIN_DIR)/grid_overlap \
     $(BIN_DIR)/cal_read_depth_from_bcd21 $(BIN_DIR)/cal_barcode_depth_from_bcd21 \
     $(BIN_DIR)/cal_twin_win_bcd_cnt $(BIN_DIR)/cal_centroid_from_read_depth \
     $(BIN_DIR)/cal_2d_overlapping_barcodes $(BIN_DIR)/small_deletion_detection \
     $(BIN_DIR)/cnv_detection

# Ensure bin/ and lib/ directories exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# Build pigz and move its binaries into bin/
pigz: | $(BIN_DIR)
	$(MAKE) -C pigz
	mv pigz/pigz $(BIN_DIR)/
	mv pigz/unpigz $(BIN_DIR)/

# Build fermikit (its own Makefile takes care of its targets)
fermikit:
	$(MAKE) -C fermikit

# Build htslib and copy its library to src/lib/
$(LIB_DIR)/libhts.a: | $(LIB_DIR)
	cd $(HTSLIB_DIR) && ./configure && $(MAKE) && $(MAKE) check
	cp $(HTSLIB_DIR)/libhts.a $(LIB_DIR)/

# Compile the C++ programs (most of these don’t need htslib so depend only when needed)

$(BIN_DIR)/cluster_reads: $(SRC_DIR)/cluster_reads.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z -pthread

$(BIN_DIR)/extract_barcode_info: $(SRC_DIR)/extract_barcode_info.cpp | $(BIN_DIR) $(LIB_DIR)/libhts.a
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $< -l hts -l z -l m -l pthread

$(BIN_DIR)/output_bam_coreinfo: $(SRC_DIR)/output_bam_coreinfo.cpp | $(BIN_DIR) $(LIB_DIR)/libhts.a
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $< -l hts -l z -l m -l pthread

$(BIN_DIR)/remove_sparse_nodes: $(SRC_DIR)/remove_sparse_nodes.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -lz

$(BIN_DIR)/cal_hap_read_depth_from_bcd21: $(SRC_DIR)/cal_hap_read_depth_from_bcd21.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/grid_overlap: $(SRC_DIR)/grid_overlap.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/cal_read_depth_from_bcd21: $(SRC_DIR)/cal_read_depth_from_bcd21.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/cal_barcode_depth_from_bcd21: $(SRC_DIR)/cal_barcode_depth_from_bcd21.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/cal_twin_win_bcd_cnt: $(SRC_DIR)/cal_twin_win_bcd_cnt.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/cal_centroid_from_read_depth: $(SRC_DIR)/cal_centroid_from_read_depth.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/cal_2d_overlapping_barcodes: $(SRC_DIR)/cal_2d_overlapping_barcodes.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/small_deletion_detection: $(SRC_DIR)/small_deletion_detection.cpp $(SRC_DIR)/cnv.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

$(BIN_DIR)/cnv_detection: $(SRC_DIR)/cnv_detection.cpp $(SRC_DIR)/cnv.cpp $(SRC_DIR)/tk.cpp $(SRC_DIR)/cgranges.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) -L $(LIB_DIR) -o $@ $^ -l z

# Clean up all built files.
clean:
	$(MAKE) -C pigz clean
	$(MAKE) -C fermikit clean
	$(MAKE) -C $(HTSLIB_DIR) clean
	rm -rf $(BIN_DIR) $(LIB_DIR)
