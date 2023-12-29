#!/bin/bash

echo "$OSTYPE"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  # Linux
	sudo apt-get update -qq
  sudo apt-get install -y --no-install-recommends libgl1-mesa-dev
	sudo apt-get install clang-format
elif [[ "$OSTYPE" == "darwin"* ]]; then
  # Mac OSX
	xcode-select --install
	brew install clang-format
else
	echo "unsupported operating system $OSTYPE"
fi
