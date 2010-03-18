################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../arith/GFE.cpp \
../arith/ModPoly.cpp \
../arith/Poly.cpp 

OBJS += \
./arith/GFE.o \
./arith/ModPoly.o \
./arith/Poly.o 

CPP_DEPS += \
./arith/GFE.d \
./arith/ModPoly.d \
./arith/Poly.d 


# Each subdirectory must supply rules for building sources it contributes
arith/%.o: ../arith/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


