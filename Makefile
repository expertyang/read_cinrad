PACKAGE = read_cinrad
default:
	cd src;make;cd ..;ln -sf src/read_cinrad .
clean:
	cd src;make clean;cd ..; rm -f read_cinrad
tar:
	tar -cvf $(PACKAGE).tar src/Makefile src/*.f90 Makefile radar_info.txt
	mv $(PACKAGE).tar $(PACKAGE)-`date +%Y%m%d%H%M`.tar
