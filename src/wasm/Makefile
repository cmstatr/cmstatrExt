

SHELL=/bin/bash

SRC=\
	webcmstatr.cpp #\
	#catch_amalgamated.cpp
	
#webcmstatr.o: $(SRC)
#	clang++ \
#		$(SRC) \
#		-std=c++17 \
#		-D WASM \
#		-o webcmstatr.o \
#		-O3 \
#		-Wcast-align \
#		-Wover-aligned \


webcmstatr.js: $(SRC)
	( \
		cd emsdk ; \
		source ./emsdk_env.sh ; \
		cd .. ; \
		emcc \
			$(SRC) \
			-std=c++17 \
			-D WASM \
			-o webcmstatr.js \
			-O3 \
			-s ASSERTIONS=1 \
			-s ALLOW_MEMORY_GROWTH=1 \
			-Wcast-align \
			-Wover-aligned \
			-sSTACK_SIZE=13107200 \
			-s "EXPORTED_RUNTIME_METHODS=['ccall']" \
			-s EXPORTED_FUNCTIONS="['_malloc', '_free', '_main', '_do_tests']" \
	)

devserver:
	python3 -m http.server

.PHONY: devserver
