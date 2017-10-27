################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/modules/BarberCushman/NutrientUptake.cpp \
../src/modules/BarberCushman/NutrientUptake2.cpp 

OBJS += \
./src/modules/BarberCushman/NutrientUptake.o \
./src/modules/BarberCushman/NutrientUptake2.o 

CPP_DEPS += \
./src/modules/BarberCushman/NutrientUptake.d \
./src/modules/BarberCushman/NutrientUptake2.d 


# Each subdirectory must supply rules for building sources it contributes
src/modules/BarberCushman/%.o: ../src/modules/BarberCushman/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


