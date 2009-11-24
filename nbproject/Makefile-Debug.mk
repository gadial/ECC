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
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_CONF=Debug
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/ecprime.o \
	${OBJECTDIR}/ellipticcurve.o \
	${OBJECTDIR}/ecbinary.o \
	${OBJECTDIR}/tests.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/coordinates.o \
	${OBJECTDIR}/curvesnisttest.o \
	${OBJECTDIR}/primes.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-lgmpxx -lgmp -lcppunit -ldl
CXXFLAGS=-lgmpxx -lgmp -lcppunit -ldl

# Fortran Compiler Flags
FFLAGS=-lgmpxx -lgmp

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Debug.mk dist/Debug/GNU-Linux-x86/ecc

dist/Debug/GNU-Linux-x86/ecc: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ecc ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/ecprime.o: nbproject/Makefile-${CND_CONF}.mk ecprime.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/ecprime.o ecprime.cpp

${OBJECTDIR}/ellipticcurve.o: nbproject/Makefile-${CND_CONF}.mk ellipticcurve.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/ellipticcurve.o ellipticcurve.cpp

${OBJECTDIR}/ecbinary.o: nbproject/Makefile-${CND_CONF}.mk ecbinary.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/ecbinary.o ecbinary.cpp

${OBJECTDIR}/tests.o: nbproject/Makefile-${CND_CONF}.mk tests.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/tests.o tests.cpp

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/coordinates.o: nbproject/Makefile-${CND_CONF}.mk coordinates.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/coordinates.o coordinates.cpp

${OBJECTDIR}/curvesnisttest.o: nbproject/Makefile-${CND_CONF}.mk curvesnisttest.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/curvesnisttest.o curvesnisttest.cpp

${OBJECTDIR}/primes.o: nbproject/Makefile-${CND_CONF}.mk primes.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/primes.o primes.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Debug
	${RM} dist/Debug/GNU-Linux-x86/ecc

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
