CC := gcc
CFLAGS := -Wall -Wextra -O2 -std=c11
CPPFLAGS := -Iinclude
LDLIBS := -lm
BIN_DIR := bin

SCALAR_STUDY := $(BIN_DIR)/scalar_study
SCALAR_GENERAL := $(BIN_DIR)/scalar_general
BLOCK_STUDY := $(BIN_DIR)/block_study
BLOCK_GENERAL := $(BIN_DIR)/block_general

.PHONY: all clean

all: $(SCALAR_STUDY) $(SCALAR_GENERAL) $(BLOCK_STUDY) $(BLOCK_GENERAL)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(SCALAR_STUDY): src/scalar/tridiag_lu_scalaire_etude.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

$(SCALAR_GENERAL): src/scalar/tridiag_lu_scalaire_general.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

$(BLOCK_STUDY): src/block/tridiag_lu_blocs_etude.c src/block/matrices_tools.c include/matrices_tools.h | $(BIN_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) src/block/tridiag_lu_blocs_etude.c src/block/matrices_tools.c -o $@ $(LDLIBS)

$(BLOCK_GENERAL): src/block/tridiag_lu_blocs_general.c src/block/matrices_tools.c include/matrices_tools.h | $(BIN_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) src/block/tridiag_lu_blocs_general.c src/block/matrices_tools.c -o $@ $(LDLIBS)

clean:
	rm -rf $(BIN_DIR)
