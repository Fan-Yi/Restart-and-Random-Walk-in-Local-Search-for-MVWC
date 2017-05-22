CC = g++
CFLAGS = -DNDEBU -O3
target = LSCC
all : $(target)

$(target) : LSCC.cpp
	$(CC) $(CFLAGS) LSCC.cpp -o $(target)

clean :
	-rm $(target)
