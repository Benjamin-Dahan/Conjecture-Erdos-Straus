CC ?= gcc
CFLAGS ?= -Wall -Wextra -Werror -std=c11 -O3

BUILD_DIR := build
TARGET := $(BUILD_DIR)/erdos_straus_sieve

.PHONY: all clean test verify

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): src/c/erdos_straus_sieve.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $< -o $@

test: $(TARGET)
	python3 tests/test_c_cli.py
	python3 src/python/verify_paper1.py $(TARGET)

verify: $(TARGET)
	./$(TARGET) --find-H 10000000
	./$(TARGET) --verify-range 500000 --depth 64

clean:
	rm -rf $(BUILD_DIR)
