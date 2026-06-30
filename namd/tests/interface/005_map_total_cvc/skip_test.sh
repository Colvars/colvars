#!/bin/bash
ARCH=$(arch)
if [ "$ARCH" = "aarch64" ]; then
  exit 0
else
  exit 1
fi
