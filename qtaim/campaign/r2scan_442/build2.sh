#!/bin/bash
cd /workspace/qe-7.5-libxc
echo "=== build start $(date) ==="
make -j 16 pw 2>&1
echo "=== build exit=$? $(date) ==="
