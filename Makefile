CC = gcc  -O0 -I ./include  -L ./lib/

all: extract_barcode sort_barcode
extract_barcode: extract_barcode.o
	$(CC) -o extract_barcode extract_barcode.o -l hts -l z -l m -l pthread

sort_barcode: sort_barcode.o utils.o
	$(CC) -o sort_barcode sort_barcode.o utils.o -I.

clean:
	rm -f *.o  

