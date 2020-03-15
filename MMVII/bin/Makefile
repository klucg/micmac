# =========== Directories
MMV1_DIR=../..
MMV1_LIBDIR=${MMV1_DIR}/lib
MMV2_DIR=${MMV1_DIR}/MMVII
MMV2_BINDIR=${MMV2_DIR}/bin
MMV2_SRCDIR=${MMV2_DIR}/src
MMV2_EXE=${MMV2_BINDIR}/MMVII
MMV2_INSTALL_PATH=${dir ${abspath ${MMV2_EXE}}}


# =========== Includes & Libraries
MMV1_INCLUDES=-I ${MMV1_DIR}/include -I${MMV1_DIR}
MMV2_INCLUDES=-I ${MMV2_DIR}/include

MMV2_LIBNAME=P2007
MMV2_LIB=${MMV2_BINDIR}/lib${MMV2_LIBNAME}.a
MMV1_LIBS=${MMV1_LIBDIR}/libelise.a ${MMV1_LIBDIR}/libANN.a
BOOST_LIBS= -lboost_system -lboost_serialization -lboost_regex -lboost_filesystem
QT_LIBS=-lXext -lGLU -lGL -lQt5Core -lQt5Gui -lQt5Xml -lQt5OpenGL -lQt5Concurrent -lQt5Widgets
EXT_LIBS=-lpthread -lX11 ${BOOST_LIBS} ${QT_LIBS}


# =========== Compiler & Flags
CXX=g++
LD=g++
CXXFLAGS= "-std=c++14" "-Wall"  "-Werror" "-O4" "-march=native" ${MMV2_INCLUDES} ${MMV1_INCLUDES}
LDFLAGS=


# ===========  Use ccache if found
CCACHE:=${shell which ccache}
ifneq (,${CCACHE})
  CXX:=${CCACHE} ${CXX}
endif


# ============ Targets
all: ${MMV2_EXE}

clean:
	rm -f ${OBJS} ${DEPS} ${MMV2_EXE} ${MMV2_LIB}


# =========== Build Lib
SRCS=${wildcard ${MMV2_SRCDIR}/*/*.cpp}
OBJS=${SRCS:.cpp=.o}
DEPS=${join ${addsuffix ., ${dir ${OBJS}}}, ${notdir ${OBJS:.o=.d}}}

%.o: %.cpp Makefile
	${CXX} -c  ${CXXFLAGS} $< -o $@

${MMV2_LIB}: ${OBJS}
	rm -f $@
	ar crs $@ $^


# ========== Build Main
MAIN_SRC=${MMV2_SRCDIR}/main.cpp
MAIN_OBJ=${MMV2_SRCDIR}/main.o


${MMV2_EXE}: ${MAIN_SRC} ${MMV2_LIB} ${MMV1_LIBS}
	$(CXX) -c $(CXXFLAGS) $(MAIN_SRC) -o $(MAIN_OBJ) -D MMVII_INSTALL_PATH='"$(MMV2_INSTALL_PATH)"'
	${CXX} ${LDFLAGS} ${MAIN_OBJ} ${MMV2_LIB} ${MMV1_LIBS} ${EXT_LIBS}  -o $@


# =========== Auto dependancies
.%.d: %.cpp
	@set -e; rm -f $@; \
	${CXX} -MM -MP -MT '${patsubst %.cpp,%.o,$<} $@' -MF $@ ${CXXFLAGS} $<   2>/dev/null

ifneq (${MAKECMDGOALS},clean)
-include ${DEPS}
endif
