dev:
	g++ -std=c++11 main.cpp system.cpp box.cpp dump.cpp particle.cpp molecule.cpp io.cpp init_position.cpp thermo.cpp boundary/*.cpp moves/*.cpp sampler/*.cpp forcefield/*.cpp rng/*.cpp -I boundary -I sampler -I moves -I forcefield -I rng -I thermo -o main.out -O3 

parser:
	g++ -std=c++11 main_parser.cpp system.cpp parser.cpp particle.cpp molecule.cpp dump.cpp box.cpp io.cpp init_position.cpp thermo.cpp boundary/*.cpp moves/*.cpp sampler/*.cpp forcefield/*.cpp rng/*.cpp -I boundary -I sampler -I moves -I forcefield -I rng -I thermo -o avbmc -O3 
	
run:
	./main.out
