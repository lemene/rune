CXX := g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -pthread
CPPFLAGS ?= -Isrc -MMD -MP

BUILD_DIR := build
BIN_DIR := $(BUILD_DIR)/bin
OBJ_DIR := $(BUILD_DIR)/obj

SRCS := src/rune.cpp
OBJS := $(patsubst src/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)
PY_SCRIPTS := $(wildcard scripts/*.py)

ifeq ($(OS),Windows_NT)
SHELL := cmd.exe
TARGET := $(BIN_DIR)/rune.exe
PY_GLOB := scripts\*.py
MKDIR = mkdir "$(subst /,\,$(1))" 2>NUL || ver > NUL
RM_RF = rmdir /s /q "$(subst /,\,$(1))" 2>NUL || ver > NUL
RM_F  = del /q "$(subst /,\,$(1))" 2>NUL || ver > NUL
COPY_PY = copy /Y $(1) $(subst /,\,$(2)/) > NUL
else
SHELL := /bin/sh
TARGET := $(BIN_DIR)/rune
PY_GLOB := scripts/*.py
MKDIR = mkdir -p $(1)
RM_RF = rm -rf $(1)
RM_F  = rm -f $(1)
COPY_PY = cp $(1) $(2)/
endif

.PHONY: all clean bin scripts

all: $(TARGET) scripts

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(OBJS) -o $@ $(CXXFLAGS)

$(OBJ_DIR)/%.o: src/%.cpp | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

scripts: $(PY_SCRIPTS) | $(BIN_DIR)
	$(call COPY_PY,$(PY_GLOB),$(BIN_DIR))

$(BIN_DIR):
	$(call MKDIR,$@)

$(OBJ_DIR):
	$(call MKDIR,$@)

clean:
	$(call RM_RF,$(BUILD_DIR))
	$(call RM_F,rune)
	$(call RM_F,rune.exe)

bin: all

-include $(DEPS)
