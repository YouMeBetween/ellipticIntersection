ellipticIntersection.exe:ellipticIntersection.cpp
ifeq ($(DEBUG), 1)
	g++ -g -o ellipticIntersection.exe ellipticIntersection.cpp
else
	g++ -o ellipticIntersection.exe ellipticIntersection.cpp
endif
