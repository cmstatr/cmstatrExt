

SHELL=/bin/bash

SRC=\
	webcmstatrExt.cpp \
	power_sim.cpp \
	../root.cpp \
	../integration.cpp \
	../acceptance.cpp

NMATH_SRC=\
	nmath/bd0.cpp nmath/chebyshev.cpp nmath/cospi.cpp \
	nmath/d1mach.cpp nmath/dnorm.cpp nmath/dpois.cpp \
	nmath/dt.cpp nmath/fmax2.cpp nmath/fmin2.cpp \
	nmath/gamma.cpp nmath/i1mach.cpp nmath/lbeta.cpp \
	nmath/lgammacor.cpp nmath/lgamma.cpp nmath/pbeta.cpp \
	nmath/pgamma.cpp nmath/pnorm.cpp nmath/pt.cpp \
	nmath/qnorm.cpp nmath/qt.cpp nmath/stirlerr.cpp \
	nmath/toms708.cpp nmath/rnorm.cpp nmath/snorm.cpp

HEADERS=\
	catch.hpp \
	testthat_catch.h \
	power_sim.h \
	../root.h \
	../integration.h \
	../acceptance.h
	
#webcmstatr.o: $(SRC)
#	clang++ \
#		$(SRC) \
#		-std=c++17 \
#		-D WASM \
#		-o webcmstatr.o \
#		-O3 \
#		-Wcast-align \
#		-Wover-aligned \


webcmstatrExt.js: $(SRC) $(HEADERS) $(NMATH_SRC)
	( \
		cd emsdk ; \
		source ./emsdk_env.sh ; \
		cd .. ; \
		emcc \
			$(SRC) \
			$(NMATH_SRC) \
			-std=c++17 \
			-D WASM \
			-o webcmstatrExt.js \
			-O3 \
			-s ASSERTIONS=1 \
			-s ALLOW_MEMORY_GROWTH=1 \
			-Wcast-align \
			-Wover-aligned \
			-sSTACK_SIZE=1310720 \
			-s "EXPORTED_RUNTIME_METHODS=['ccall']" \
			-s EXPORTED_FUNCTIONS="['_malloc', '_free', '_main', '_do_tests', \
			                        '_k_equiv_two_sample', '_power_mean', \
			                        '_p_equiv_two_sample']" \
	)

devserver:
	python3 -m http.server

.PHONY: devserver
