PACKAGE = read_cinrad
default:
	cd src;make;cd ..;ln -sf src/read_cinrad .;ln -sf src/proc_radar .
clean:
	rm -f *~;make dataclean;cd src;make clean;cd ..; rm -f read_cinrad proc_radar
dataclean:
	rm -f fort.* raw* qced* limit* remap* radar_info_in.txt radar_info_out.txt radar*dist.txt radar*used.txt *.3dv
tar:
	tar -cvf $(PACKAGE).tar src/Makefile src/*.f90 Makefile radar_info.txt namelist.input
	mv $(PACKAGE).tar $(PACKAGE)-`date +%Y%m%d%H%M`.tar
