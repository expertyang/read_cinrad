PACKAGE = read_cinrad
default:
	cd src;make;cd ..;ln -sf src/read_cinrad .;ln -sf src/proc_radar .
clean:
	rm -f raw* *~;cd src;make clean; rm -f *~;cd ..; rm -f read_cinrad proc_radar
tar:
	tar -cvf $(PACKAGE).tar src/Makefile src/*.f90 Makefile radar_info.txt namelist.input
	mv $(PACKAGE).tar $(PACKAGE)-`date +%Y%m%d%H%M`.tar
