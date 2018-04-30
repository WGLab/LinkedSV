CC = gcc   -O0 -g -I ./include -L ./lib/

SOBJ = extract_barcode.o 

EXTRACT_BARCODE: $(SOBJ)
	$(CC)    -o extract_barcode $(SOBJ) -l hts -l z -l m -l pthread
clean:
	rm -f *.o  
