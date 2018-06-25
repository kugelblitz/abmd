#!/bin/bash
/Users/dan/Soft/valgrind/bin/valgrind --track-origins=yes ./cmake-build-debug/dde   2>&1 | less
