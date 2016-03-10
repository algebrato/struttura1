LIBS	:=-lblas -llapack -llapacke
FLAGS	:=-O3
CC	:= gcc
OBJ	:= Bands.o
EXE	:= pw.x


%.o : %.c
	$(CC) ${LIBS} ${FLAGS}  -c $< -o $@

${EXE} : ${OBJ}
	$(CC) ${LIBS} ${FLAGS} -o $@ $^

.PHONY: clean
clean:
	rm -rf *.o *.x *.out *.dat


