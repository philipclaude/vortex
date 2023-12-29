#!/bin/sh

clang-format --dry-run -Werror src/*.cpp src/*.h src/*.hpp src/math/*.h src/math/*.hpp src/math/*.cpp
