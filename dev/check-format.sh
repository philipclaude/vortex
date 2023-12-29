#!/bin/sh

clang-format --verbose --style=Google --dry-run -Werror src/*.cpp src/*.h src/*.hpp src/math/*.h src/math/*.hpp src/math/*.cpp
