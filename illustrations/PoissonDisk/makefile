include ../../Makefile.in

main: main.o disk.o gui.o gui.moc.o
	$(CXX) $(CXXFLAGS) main.o disk.o gui.o gui.moc.o -F/usr/local/lib -framework QtCore -framework QtGui -L$(DIR_LIB) $(LIB) -o $@ 

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -I/usr/local/include/QtCore -I/usr/local/include/QtGui -c $< -o $@

disk.o: PoissonDisk.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -c PoissonDisk.cpp -o $@

gui.o: gui.cpp
	$(CXX) $(CXXFLAGS) $(DIR_INC) -I/usr/local/include/QtCore -I/usr/local/include/QtGui -c gui.cpp -o $@

gui.moc.o: gui.moc.cpp
	g++ -std-c++11 -c -I/usr/local/include/QtCore -I/usr/local/include/QtGui gui.moc.cpp -o $@

%.moc.cpp: %.h
		moc $< -o $@

.PHONY:clean
clean:
	rm -f *.o