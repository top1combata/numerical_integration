BUILD_DIR   := bin
COMPILER    := g++
SOURCES     := main.cpp functions.cpp utils.cpp
FLAGS       := -Wall -O2 -std=c++20
DEBUG_FLAGS := -O0 -g -std=c++20 
EXEC_NAME   := main
LIBS        := -lgmp


all: build_exec

debug: build_dir
	@${COMPILER} ${SOURCES} ${DEBUG_FLAGS} ${LIBS} -o ${BUILD_DIR}/debug

build_dir:
	@mkdir -p ${BUILD_DIR}

build_exec: build_dir
	@${COMPILER} ${SOURCES} ${FLAGS} ${LIBS} -o ${BUILD_DIR}/${EXEC_NAME}

run: build_exec
	@./${BUILD_DIR}/main

clean:
	@rm -rf ${BUILD_DIR}
