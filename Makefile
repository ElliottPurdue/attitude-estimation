# Host build for the test suite and simulation driver.
#
# -Wall -Wextra -Werror because this library is meant to be dropped into
# firmware, where a warning is usually a defect that has not surfaced yet.
# C99 for the toolchain floor; nothing here needs anything newer.

# GNU Make predefines CC as "cc", which counts as set, so ?= would leave it
# alone and MinGW has no cc. Replace only make's own default, so CC=clang on the
# command line still wins.
ifeq ($(origin CC),default)
CC := gcc
endif

CFLAGS  := -std=c99 -O2 -Wall -Wextra -Werror -pedantic
LDFLAGS := -lm

SRC     := $(wildcard src/*.c)
TESTS   := $(wildcard tests/*.c)

BUILD   := build
TEST_BIN  := $(BUILD)/run_tests
BENCH_BIN := $(BUILD)/run_filter
TIME_BIN  := $(BUILD)/bench_cost

.PHONY: all test bench time clean

all: test

$(BUILD):
	@mkdir -p $(BUILD)

$(TEST_BIN): $(SRC) $(TESTS) | $(BUILD)
	$(CC) $(CFLAGS) $(SRC) $(TESTS) -o $@ $(LDFLAGS)

test: $(TEST_BIN)
	@./$(TEST_BIN)

$(BENCH_BIN): $(SRC) bench/run_filter.c | $(BUILD)
	$(CC) $(CFLAGS) $(SRC) bench/run_filter.c -o $@ $(LDFLAGS)

bench: $(BENCH_BIN)

$(TIME_BIN): $(SRC) bench/bench_cost.c | $(BUILD)
	$(CC) $(CFLAGS) $(SRC) bench/bench_cost.c -o $@ $(LDFLAGS)

time: $(TIME_BIN)
	@./$(TIME_BIN)

clean:
	@rm -rf $(BUILD)
