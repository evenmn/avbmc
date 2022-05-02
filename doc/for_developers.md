# Build with g++ (for developers)
If you are a developer, you might want to build the code using g++ and run various debugging tools. To build release, run 
``` bash
make release
```
It will build the system specified in `src/main.cpp`, which is never used when building from source using `setuptools` (`pip install avbmc`). To use Valgrind for memory checks, first build the code with the debugging tools,
``` bash
make debug
```
and then run the code with `Valgrind`,
``` bash
make valgrind
```
and inspect `valgrind-out.txt`. This works on Linux only. For a general code check, run
``` bash
make cppcheck
```
