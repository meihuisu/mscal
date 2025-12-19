# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} 
AM_LDFLAGS = ${LDFLAGS} -lm


TARGETS = libmscal.a libmscal.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libmscal.so ${prefix}/lib
	cp libmscal.a ${prefix}/lib
	cp mscal.h ${prefix}/include

clean:
	rm -rf $(TARGETS)
	rm -rf *.o

libmscal.a: mscal_static.o
	$(AR) rcs $@ $^

libmscal.so: mscal.o 
	$(CC) -shared $(AM_FCFLAGS) -o libmscal.so $^ $(AM_LDFLAGS)

mscal.o: mscal.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

mscal_static.o: mscal.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

