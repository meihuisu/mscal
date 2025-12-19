# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
if WITH_MSCAL_LARGEDATA_DIR
	echo "USING WITH_MSCAL_LARGEDATA_DIR in mscal"
else
	echo "NOT USING WITH_MSCAL_LARGEDATA_DIR in mscal"
	if [ ! -d "mscal" ]; then \
		./make_data_files.py; \
	fi 
endif

install:
if WITH_MSCAL_LARGEDATA_DIR
	if [ ! -d "${prefix}/data/mscal" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		$(LN_S) -f $(MSCAL_LARGEDATA_DIR) ${prefix}/data/mscal; \
		echo "USING WITH_MSCAL_LARGEDATA_DIR in mscal"; \
	fi
	cp config ${prefix}/data
else
	if [ ! -d "${prefix}/data/mscal" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		mv mscal ${prefix}/data; \
		$(LN_S) -f ${prefix}/data/mscal .; \
		echo "NOT USING WITH_MSCAL_LARGEDATA_DIR in mscal"; \
	fi
	cp config ${prefix}/data
endif

clean:
	echo "clean"
