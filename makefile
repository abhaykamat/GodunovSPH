grpscheme: main.cpp mesh.o flow.o
	g++ -o grpscheme main.cpp mesh.o

mesh.o: mesh/mesh.cpp flow.o
	g++ -c mesh/mesh.cpp flow.o

flow.o: flow/flow.cpp
	g++ -c flow/flow.cpp