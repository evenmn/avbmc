dev:
	mpicxx -std=c++11 src/main.cpp src/system.cpp src/dump.cpp src/box.cpp src/particle.cpp src/io.cpp src/init_position.cpp src/thermo.cpp src/boundary/*.cpp src/moves/*.cpp src/sampler/*.cpp src/forcefield/*.cpp src/rng/*.cpp -I src/ -I src/boundary -I src/sampler -I src/moves -I src/forcefield -I src/rng -I src/thermo -o main.out -O3

parser:
	mpicxx -std=c++11 src/main_parser.cpp src/system.cpp src/parser.cpp src/particle.cpp src/dump.cpp src/box.cpp src/io.cpp src/init_position.cpp src/thermo.cpp src/boundary/*.cpp src/moves/*.cpp src/sampler/*.cpp src/forcefield/*.cpp src/rng/*.cpp -I src -I src/boundary -I src/sampler -I src/moves -I src/forcefield -I src/rng -I src/thermo -o avbmc -O3 
	
debug:
	mpicxx -std=c++11 -Wall -Wextra -ggdb3 src/main3.cpp src/system.cpp src/dump.cpp src/box.cpp src/particle.cpp src/io.cpp src/init_position.cpp src/thermo.cpp src/boundary/*.cpp src/moves/*.cpp src/sampler/*.cpp src/forcefield/*.cpp src/rng/*.cpp -I src/ -I src/boundary -I src/sampler -I src/moves -I src/forcefield -I src/rng -I src/thermo -o main.out -O3

run:
	./main.out

valgrind:
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./main.out
