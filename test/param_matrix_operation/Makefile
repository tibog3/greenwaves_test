# User Test
#------------------------------------
APP              = test
# App sources
APP_SRCS         = matrix_operation.c
# App includes
APP_INC	         =
# Linker flags
APP_LDFLAGS      =
# Custom linker
APP_LINK_SCRIPT  =


MAT_SIZE ?= 64

# Utilisez la variable MAT_SIZE dans APP_CFLAGS
APP_CFLAGS += -DMAT_SIZE=$(MAT_SIZE)

include $(RULES_DIR)/pmsis_rules.mk
