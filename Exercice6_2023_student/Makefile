CC = g++
CFLAGS = -std=c++11 -Wall -pedantic -g
EXEC_NAME = Exercice6_2023_solution
EXEC_NAME_STD = Exercice6_2023_student
INCLUDES =
LIBS =
OBJ_FILES = Exercice6_2023_solution.o
OBJ_FILES_STD = Exercice6_2023_student.o
EXEC_SUFFIX = .exe

# add suffix to exec
EXEC_NAME_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME))
EXEC_NAME_STD_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME_STD))

all : solution student

clean :
	rm $(EXEC_NAME_SUFFIX) $(EXEC_NAME_STD_SUFFIX)  $(OBJ_FILES) *.out
 
solution : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME_SUFFIX) $(OBJ_FILES) $(LIBS)

student : $(OBJ_FILES_STD)
	$(CC) -o $(EXEC_NAME_STD_SUFFIX) $(OBJ_FILES_STD) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

