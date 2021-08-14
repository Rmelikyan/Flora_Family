export OPENGL=0

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../ -Wl,-rpath,./ problem.c -L. -lrebound -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound:
	@-rm -f librebound.so
	@ln -s ../librebound.so .

clean:
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
