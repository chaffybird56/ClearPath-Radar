################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Each subdirectory must supply rules for building sources it contributes
RadarProcessor_lab4.obj: ../RadarProcessor_lab4.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C6000 Compiler'
	"C:/ti/ccsv5/tools/compiler/c6000_7.4.4/bin/cl6x" --abi=coffabi -g --include_path="C:/ti/ccsv5/tools/compiler/c6000_7.4.4/include" --include_path="Y:/COE4TL4/dsplib/include" --include_path="Y:/COE4TL4/C6xCSL/include" --include_path="Y:/DSK6713/c6000/dsk6713/include" --include_path="Y:/COE4TL4/RadarProcessor/include" --include_path="Y:/COE4TL4/RadarProcessor/data/group15/no_barrier" --display_error_number --diag_warning=225 --diag_wrap=off --mem_model:const=far --mem_model:data=far --preproc_with_compile --preproc_dependency="RadarProcessor_lab4.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

Vectors_intr.obj: ../Vectors_intr.asm $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C6000 Compiler'
	"C:/ti/ccsv5/tools/compiler/c6000_7.4.4/bin/cl6x" --abi=coffabi -g --include_path="C:/ti/ccsv5/tools/compiler/c6000_7.4.4/include" --include_path="Y:/COE4TL4/dsplib/include" --include_path="Y:/COE4TL4/C6xCSL/include" --include_path="Y:/DSK6713/c6000/dsk6713/include" --include_path="Y:/COE4TL4/RadarProcessor/include" --include_path="Y:/COE4TL4/RadarProcessor/data/group15/no_barrier" --display_error_number --diag_warning=225 --diag_wrap=off --mem_model:const=far --mem_model:data=far --preproc_with_compile --preproc_dependency="Vectors_intr.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '


