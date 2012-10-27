grpscheme: main.cpp scheme.o riemann.o mesh.o flow.o 
	clang++ -o grpscheme main.cpp scheme.o riemann.o mesh.o flow.o -arch x86_64

scheme.o: scheme/scheme.cpp
	clang++ -c scheme/scheme.cpp -arch x86_64

riemann.o: riemann/riemann.cpp
	clang++ -c riemann/riemann.cpp -arch x86_64

mesh.o: mesh/mesh.cpp
	clang++ -c mesh/mesh.cpp -arch x86_64

flow.o: flow/flow.cpp
	clang++ -c flow/flow.cpp -arch x86_64