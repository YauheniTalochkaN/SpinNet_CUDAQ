.PHONY: all clean

NVQPP = nvq++

NVQPPFLAGS = -std=c++20 -Wall -O2 -fPIC --target dynamics

PROJECT_DIR = .
CUDAQ_LIB_DIR = /opt/nvidia/cudaq/lib
CUDAQ_INCLUDE_DIR = /opt/nvidia/cudaq/include

LIBS = -L$(CUDAQ_LIB_DIR) -lcudaq -Wl,-rpath,$(CUDAQ_LIB_DIR)

INCLUDES = -I$(PROJECT_DIR)/include -I$(CUDAQ_INCLUDE_DIR)

SRCS = $(shell find $(PROJECT_DIR)/src/ -name '*.cc') \
       $(PROJECT_DIR)/SpinDynamics.cc
	   
OBJS = $(patsubst %.cc,%.o,$(SRCS))

GREEN := $(shell tput -Txterm setaf 2)
BLUE  := $(shell tput -Txterm setaf 6)
RESET := $(shell tput -Txterm sgr0)

all: SpinDynamics

clean:
	@rm -rf SpinDynamics $(OBJS)

SpinDynamics: $(OBJS)
	$(info $(GREEN)Linking all objects and libs into SpinDynamics$(RESET))
	@$(NVQPP) $(NVQPPFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)
	@rm -rf $(OBJS)
	$(info $(BLUE)Built the target SpinDynamics$(RESET))

%.o: %.cc
	$(info $(GREEN)Compiling $< into $@ $(RESET))
	@$(NVQPP) $(NVQPPFLAGS) $(INCLUDES) -c $< -o $@