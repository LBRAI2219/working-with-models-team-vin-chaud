################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/export/Text/TabledOutput.cpp \
../src/export/Text/WriteModel2File.cpp 

OBJS += \
./src/export/Text/TabledOutput.o \
./src/export/Text/WriteModel2File.o 

CPP_DEPS += \
./src/export/Text/TabledOutput.d \
./src/export/Text/WriteModel2File.d 


# Each subdirectory must supply rules for building sources it contributes
src/export/Text/%.o: ../src/export/Text/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


