JO_COMPILE_WITH_VIDEO_MODULE = 1
JO_COMPILE_WITH_BACKUP_MODULE = 1
JO_COMPILE_WITH_TGA_MODULE = 1
JO_COMPILE_WITH_AUDIO_MODULE = 0
JO_COMPILE_WITH_3D_MODULE = 1
JO_COMPILE_WITH_PSEUDO_MODE7_MODULE = 1
JO_COMPILE_WITH_EFFECTS_MODULE = 1
JO_COMPILE_WITH_FAST_BUT_LESS_ACCURATE_MATH = 1
JO_COMPILE_WITH_PRINTF_SUPPORT = 1
JO_DEBUG = 0
JO_NTSC = 1
JO_GLOBAL_MEMORY_SIZE_FOR_MALLOC = 65536
SRCS=main.c ssv.c newcol.c
JO_ENGINE_SRC_DIR=../../jo_engine
COMPILER_DIR=../../Compiler
include $(COMPILER_DIR)/COMMON/jo_engine_makefile
