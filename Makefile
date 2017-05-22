SHELL =         /bin/sh
CXX =		g++
MAKE =		make
########## OpenMP
#DEBUG =		-g -Wall -Wno-deprecated -DMULTICPU -fopenmp
#DEBUG =		-O4 -Wno-deprecated -DMULTICPU -fopenmp
########## Normal
#DEBUG =		-g -Wall -Wno-deprecated -frounding-math
#DEBUG =		-O2 -Wno-deprecated -frounding-math
#DEBUG =		-O3 -Wno-deprecated -frounding-math
DEBUG =		-O4 -Wno-deprecated -frounding-math
##########
#CFLAG =		$(DEBUG) -m64
CFLAG =		$(DEBUG)
##########
MACDEFS =	-D__MAC__
LINUXDEFS =	-D__LINUX__
##########
MACINC =	-I/usr/X11R6/include -I/usr/X11R6/include/freetype2
LINUXINC =	-I/usr/X11R6/include
##########
MACLIB =	`pkg-config --libs opencv` \
       		-lCGAL \
		-framework OpenGL -framework GLUT \
                -lgsl -lgslcblas -lm
LINUXLIB =	`pkg-config --libs opencv` \
       		-lCGAL \
	  	-lglut -lGLU \
                -lgsl -lgslcblas -lm
##########
SRC =           main.cpp \
		dual.cpp feature.cpp matching.cpp \
		craft.cpp project.cpp stitch.cpp search.cpp \
		mst.cpp \
		random.cpp gene.cpp \
		meshui.cpp sheetui.cpp render.cpp \
		timer.cpp
HDR =		common.h ui.h random.h gene.h timer.h
OBJ =		$(SRC:.cpp=.o)
BIN =		unfold
#
.cpp.o:	$(SRC) $(HDR)
	$(CXX) $(DEFS) $(CFLAG) $(INC) -c $<

default:
	make $(BIN)

mac	:
	make $(OBJ) "DEFS=$(MACDEFS)" "INC=$(MACINC)"
	$(CXX) -o $(BIN) $(MACDEFS) $(OBJ) $(MACLIB)

linux	:
	make $(OBJ) "DEFS=$(LINUXDEFS)" "INC=$(LINUXINC)"
	$(CXX) -o $(BIN) $(LINUXDEFS) $(OBJ) $(LINUXLIB)

clean:
	rm -f *.o *~* core $(BIN)

#
$(OBJ): $(HDR)

