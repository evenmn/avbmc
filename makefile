dev:
	g++ -std=c++11 src/main.cpp src/dump.cpp src/box.cpp src/particle.cpp src/io.cpp src/init_position.cpp src/thermo.cpp src/boundary/*.cpp src/moves/*.cpp src/sampler/*.cpp src/forcefield/*.cpp src/rng/*.cpp -I src/ -I src/boundary -I src/sampler -I src/moves -I src/forcefield -I src/rng -I src/thermo -o main.out -O3

parser:
	g++ -std=c++11 src/main_parser.cpp src/parser.cpp src/particle.cpp src/dump.cpp src/box.cpp src/io.cpp src/init_position.cpp src/thermo.cpp src/boundary/*.cpp src/moves/*.cpp src/sampler/*.cpp src/forcefield/*.cpp src/rng/*.cpp -I src -I src/boundary -I src/sampler -I src/moves -I src/forcefield -I src/rng -I src/thermo -o avbmc -O3 
	
run:
	./main.out
