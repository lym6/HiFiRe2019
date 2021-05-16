CC = /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgc++ 
CFLAGS = -c
LIBB = -L/usr/local/lib -lgeos
LIBRA = lib/libshp.a lib/librtree.a
DEBUG =

all: cpu_join.o openacc.o 
	$(CC) -fast -Minfo -Mbounds -acc -ta=tesla -mp -o spatialjoin cpu_join.o openacc.o $(LIBB) $(LIBRA)
			
cpu_join.o: cpu_join.cpp
	$(CC) $(CFLAGS) -std=c++11 -mp cpu_join.cpp $(LIBRA) $(LIBB)

openacc.o: openacc.cpp
	$(CC) $(CFLAGS) -mp -Mlarge_arrays -fast -Minfo messages -Mbounds -acc -ta=tesla openacc.cpp $(LIBB) $(LIBRA)