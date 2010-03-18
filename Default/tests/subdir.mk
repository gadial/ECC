################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tests/bincurvetest.cpp \
../tests/curvesnisttest.cpp \
../tests/padictest.cpp 

OBJS += \
./tests/bincurvetest.o \
./tests/curvesnisttest.o \
./tests/padictest.o 

CPP_DEPS += \
./tests/bincurvetest.d \
./tests/curvesnisttest.d \
./tests/padictest.d 


# Each subdirectory must supply rules for building sources it contributes
tests/%.o: ../tests/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


