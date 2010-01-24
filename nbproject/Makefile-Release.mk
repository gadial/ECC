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
CND_CONF=Release
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/ecprime.o \
	${OBJECTDIR}/adicops.o \
	${OBJECTDIR}/tests/padictest.o \
	${OBJECTDIR}/tests/curvesnisttest.o \
	${OBJECTDIR}/tests.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/coordinates.o \
	${OBJECTDIR}/primes.o \
	${OBJECTDIR}/ellipticcurve.o \
	${OBJECTDIR}/arith/Poly.o \
	${OBJECTDIR}/hcp.o \
	${OBJECTDIR}/ecbinary.o \
	${OBJECTDIR}/zp_int.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/GNU-Linux-x86/ecc

dist/Release/GNU-Linux-x86/ecc: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ecc ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/ecprime.o: nbproject/Makefile-${CND_CONF}.mk ecprime.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/ecprime.o ecprime.cpp

${OBJECTDIR}/adicops.o: nbproject/Makefile-${CND_CONF}.mk adicops.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/adicops.o adicops.cpp

${OBJECTDIR}/tests/padictest.o: nbproject/Makefile-${CND_CONF}.mk tests/padictest.cpp 
	${MKDIR} -p ${OBJECTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/tests/padictest.o tests/padictest.cpp

${OBJECTDIR}/tests/curvesnisttest.o: nbproject/Makefile-${CND_CONF}.mk tests/curvesnisttest.cpp 
	${MKDIR} -p ${OBJECTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/tests/curvesnisttest.o tests/curvesnisttest.cpp

${OBJECTDIR}/tests.o: nbproject/Makefile-${CND_CONF}.mk tests.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/tests.o tests.cpp

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/coordinates.o: nbproject/Makefile-${CND_CONF}.mk coordinates.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/coordinates.o coordinates.cpp

${OBJECTDIR}/primes.o: nbproject/Makefile-${CND_CONF}.mk primes.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/primes.o primes.cpp

${OBJECTDIR}/ellipticcurve.o: nbproject/Makefile-${CND_CONF}.mk ellipticcurve.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/ellipticcurve.o ellipticcurve.cpp

${OBJECTDIR}/arith/Poly.o: nbproject/Makefile-${CND_CONF}.mk arith/Poly.cpp 
	${MKDIR} -p ${OBJECTDIR}/arith
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/arith/Poly.o arith/Poly.cpp

${OBJECTDIR}/hcp.o: nbproject/Makefile-${CND_CONF}.mk hcp.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/hcp.o hcp.cpp

${OBJECTDIR}/ecbinary.o: nbproject/Makefile-${CND_CONF}.mk ecbinary.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/ecbinary.o ecbinary.cpp

${OBJECTDIR}/zp_int.o: nbproject/Makefile-${CND_CONF}.mk zp_int.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/zp_int.o zp_int.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/ecc

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
