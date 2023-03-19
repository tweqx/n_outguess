EXEC=n_outguess

$(EXEC): md5.o arc.o iterator.o main.o
	gcc $^ -Ljpeg-6b-steg -l:libjpeg.a -o $@

%.o: %.c
	gcc -c $< -o $@

.PHONY: clean mrproper
clean:
	rm *.o
mrproper: clean
	rm $(EXEC)
