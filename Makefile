CC = g++ 
CFLAGS = -g#-O3
LFLAGS = #
OBJECTS = main.o routines.o pimc.o path.o

#output.txt: main.exe
#	./main.exe > output.txt

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main.exe

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< 
