CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c11
LFLAGS = -lcheck -lsubunit -lm -lrt -lpthread
GCOVFLAGS = -fprofile-arcs -ftest-coverage

SRC = s21_matrix.c
TEST_SRC = s21_matrix_test.c
OBJ = $(SRC:.c=.o)
TEST_OBJ = $(TEST_SRC:.c=.o)
TARGET = s21_matrix

.PHONY: all clean test gcov_report count_depth cppcheck valgrind

all: $(TARGET).a

$(TARGET).a: $(OBJ)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(OBJ) $(TEST_OBJ)
	$(CC) $(CFLAGS) $^ -o test_runner $(LFLAGS)
	./test_runner

gcov_report: test
	$(CC) $(CFLAGS) $(GCOVFLAGS) $(SRC) $(TEST_SRC) -o gcov_runner $(LFLAGS)
	./gcov_runner
	lcov -t "gcov_report" -o gcov_report.info -c -d .
	lcov -r gcov_report.info -o gcov_report.info
	genhtml -o report gcov_report.info
#	open ./report/index.html

valgrind: test
	valgrind ./test_runner

cppcheck:
	cppcheck *.c

clean:
	rm -rf *.o *.a *.gcno *.gcda *.info test_runner gcov_runner report tests/*.o

count_depth:
	find . -name '*.c' -exec awk '\
	  { for (i=1; i<=length($$0); i++) {\
	      c = substr($$0,i,1);\
	      if (c == "{") depth++;\
	      else if (c == "}") depth--;\
	    }\
	    if (depth > max) max = depth;\
	  }\
	  END { print FILENAME ": max depth = " max }\
	' {} \;
