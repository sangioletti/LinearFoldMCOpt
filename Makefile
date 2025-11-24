################################
# Makefile
#
# author: He Zhang
# edited by: 03/2019
################################

CC=g++
DEPS=src/bpp.cpp src/LinearPartition.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h
CFLAGS=-std=c++11 -O3
ARCH_FLAGS ?=
ifeq ($(shell uname -s),Darwin)
	ARCH_FLAGS += -arch arm64 -arch x86_64
endif
UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	SHARED_EXT=dylib
	SHARED_FLAG=-dynamiclib
else
	SHARED_EXT=so
	SHARED_FLAG=-shared
endif
LIB_TARGET=liblinearpartition_v.$(SHARED_EXT)
BUILD_DIR=build
.PHONY : clean linearpartition liblinearpartition
objects=bin/linearpartition_v bin/linearpartition_c

linearpartition: src/LinearPartition.cpp $(DEPS) 
		chmod +x linearpartition draw_bpp_plot draw_heatmap
		mkdir -p bin
		$(CC) $(ARCH_FLAGS) src/LinearPartition.cpp $(CFLAGS) -Dlpv -o bin/linearpartition_v 
		$(CC) $(ARCH_FLAGS) src/LinearPartition.cpp $(CFLAGS) -o bin/linearpartition_c

liblinearpartition: $(LIB_TARGET)

$(BUILD_DIR):
		mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/LinearPartition_v.o: src/LinearPartition.cpp $(DEPS) | $(BUILD_DIR)
		$(CC) $(ARCH_FLAGS) $(CFLAGS) -fPIC -Dlpv -DLP_DISABLE_MAIN -c $< -o $@

$(BUILD_DIR)/LinearPartitionAPI_v.o: src/LinearPartitionAPI.cpp src/LinearPartitionAPI.h src/LinearPartition.h | $(BUILD_DIR)
		$(CC) $(ARCH_FLAGS) $(CFLAGS) -fPIC -Dlpv -c $< -o $@

$(LIB_TARGET): $(BUILD_DIR)/LinearPartition_v.o $(BUILD_DIR)/LinearPartitionAPI_v.o
		$(CC) $(ARCH_FLAGS) $(SHARED_FLAG) $^ -o $@

clean:
	-rm -f $(objects) $(LIB_TARGET)
	-rm -rf $(BUILD_DIR)