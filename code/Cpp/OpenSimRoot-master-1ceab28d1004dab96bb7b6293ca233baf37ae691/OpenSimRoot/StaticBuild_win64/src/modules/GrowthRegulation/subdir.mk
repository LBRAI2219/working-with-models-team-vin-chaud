################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/modules/GrowthRegulation/StressAndPlasticity.cpp 

OBJS += \
./src/modules/GrowthRegulation/StressAndPlasticity.o 

CPP_DEPS += \
./src/modules/GrowthRegulation/StressAndPlasticity.d 


# Each subdirectory must supply rules for building sources it contributes
src/modules/GrowthRegulation/%.o: ../src/modules/GrowthRegulation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	x86_64-w64-mingw32-g++ -std=c++14 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


