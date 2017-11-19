WORKER_CPP=./reference/rank.cpp
WORKER_EXE=./reference/UV/rank
TASKS_DIR=./task
MATRIX_DIR=${TASKS_DIR}/UV
SETUP_PY=${TASKS_DIR}/setup_tasks.py
TASKS_FILE=${TASKS_DIR}/tasks.txt
TEST=./your/rank

default: setup
	@

setup: compile
	@mkdir -p ${TASKS_DIR}
	@mkdir -p ${MATRIX_DIR}
	@echo Creating ${TASKS_FILE}
	@python3 ${SETUP_PY} ${MATRIX_DIR} ${WORKER_EXE} > ${TASKS_FILE}
	@echo Now run: make check

compile: 
	@echo Compiling reference worker ${WORKER_CPP}
	@g++ -std=c++11 -O3 ${WORKER_CPP} -o ${WORKER_EXE}

check:
	@[ -e ${TASKS_FILE} ] || (echo File ${TASKS_FILE} not found. Run: make setup; exit 1)
	@[ -e ${TEST} ] || (echo Executable ${TEST} not found. Compile your implementation of rank.cpp there.; exit 1)
	@while read T; do \
		f=$$(echo $$T | cut -f 1 -d ' '); \
		R=$$(echo $$T | cut -f 2 -d ' '); \
		echo ; \
		echo Testfile: $$f; \
		echo Expected: $$R; \
		r=$$(echo $$f | ${TEST} 2> /dev/null); \
		echo Received: $$r; \
	done < ${TASKS_FILE}
