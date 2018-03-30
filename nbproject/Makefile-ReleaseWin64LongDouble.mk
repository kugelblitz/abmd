#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=x86_64-w64-mingw32-gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=ReleaseWin64LongDouble
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/abm.o \
	${OBJECTDIR}/src/api.o \
	${OBJECTDIR}/src/coeffs.o \
	${OBJECTDIR}/src/poly.o \
	${OBJECTDIR}/src/queue.o


# C Compiler Flags
CFLAGS=-m64

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-static-libgcc

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dde-longdouble.dll

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dde-longdouble.dll: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	x86_64-w64-mingw32-gcc -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/dde-longdouble.dll ${OBJECTFILES} ${LDLIBSOPTIONS} -shared -fPIC

${OBJECTDIR}/src/abm.o: src/abm.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -Wall -s -DUSE_LONG_DOUBLE -I. -Isrc -std=c99 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/abm.o src/abm.c

${OBJECTDIR}/src/api.o: src/api.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -Wall -s -DUSE_LONG_DOUBLE -I. -Isrc -std=c99 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/api.o src/api.c

${OBJECTDIR}/src/coeffs.o: src/coeffs.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -Wall -s -DUSE_LONG_DOUBLE -I. -Isrc -std=c99 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/coeffs.o src/coeffs.c

${OBJECTDIR}/src/poly.o: src/poly.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -Wall -s -DUSE_LONG_DOUBLE -I. -Isrc -std=c99 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/poly.o src/poly.c

${OBJECTDIR}/src/queue.o: src/queue.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -Wall -s -DUSE_LONG_DOUBLE -I. -Isrc -std=c99 -fPIC  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/queue.o src/queue.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
