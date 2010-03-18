################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../adicops.cpp \
../challange_crack.cpp \
../cmd.cpp \
../coordinates.cpp \
../ecbinary.cpp \
../ecprime.cpp \
../elgamal.cpp \
../ellipticcurve.cpp \
../hcp.cpp \
../main.cpp \
../primes.cpp \
../tests.cpp \
../zp_int.cpp 

OBJS += \
./adicops.o \
./challange_crack.o \
./cmd.o \
./coordinates.o \
./ecbinary.o \
./ecprime.o \
./elgamal.o \
./ellipticcurve.o \
./hcp.o \
./main.o \
./primes.o \
./tests.o \
./zp_int.o 

CPP_DEPS += \
./adicops.d \
./challange_crack.d \
./cmd.d \
./coordinates.d \
./ecbinary.d \
./ecprime.d \
./elgamal.d \
./ellipticcurve.d \
./hcp.d \
./main.d \
./primes.d \
./tests.d \
./zp_int.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


