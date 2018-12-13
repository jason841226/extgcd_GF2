CC=    gcc
CXX=   g++
LD=    gcc


ifndef PROJ_DIR
#PROJ_DIR = ../ref
PROJ_DIR = ../avx2
endif

SRCS_DIR = ../srcs


CFLAGS= -O3 -std=c99 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer
CXXFLAGS= -O3 -Wall -Wextra -fno-exceptions -fno-rtti -nostdinc++
INCPATH= -I/usr/local/include -I/opt/local/include -I/usr/include -I$(SRCS_DIR) -I$(PROJ_DIR)
LDFLAGS= -fsanitize=address
LIBPATH= -L/usr/local/lib -L/opt/local/lib -L/usr/lib
LIBS=    -lcrypto -lasan


ifeq ($(PROJ_DIR),../avx2)
CFLAGS += -mavx2
CXXFLAGS += -mavx2
endif


C_SRCS = $(wildcard $(SRCS_DIR)/*.c)
C_SRCS_O = $(C_SRCS:.c=.o)
C_SRCS_O_ND = $(subst $(SRCS_DIR)/,,$(C_SRCS_O))

ARCH_SRCS = $(wildcard $(PROJ_DIR)/*.c)
ARCH_SRCS_O = $(ARCH_SRCS:.c=.o)
ARCH_SRCS_O_ND = $(subst $(PROJ_DIR)/,,$(ARCH_SRCS_O))


OBJ = $(C_SRCS_O_ND) $(ARCH_SRCS_O_ND)

EXE= extgcd-test div1step-test transmat-test gf4591-test

ifeq ($(PROJ_DIR),../avx2)
EXE += gf4591_avx2-test
endif


ifdef DEBUG
        CFLAGS+=  -D_DEBUG_ -g
        CXXFLAGS+= -D_DEBUG_ -g
endif

ifdef NO_SSE
	CFLAGS += -D_NO_SSE_
	CXXFLAGS += -D_NO_SSE_
endif

ifdef AVX2
	CFLAGS += -mavx2 -D_USE_AVX2_
	CXXFLAGS += -mavx2 -D_USE_AVX2_
endif

ifdef AVX
	CFLAGS += -mavx -D_USE_AVX_
	CXXFLAGS += -mavx -D_USE_AVX_
endif

ifdef GPROF
	CFLAGS += -pg
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

.PHONY: all tests tables clean

all: $(OBJ) $(EXE)


%-test: $(OBJ) %-test.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%-benchmark: $(OBJ) %-benchmark.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

%.o: $(PROJ_DIR)/%.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

%.o: $(SRCS_DIR)/%.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $<

tests:
	cd tests; make

tables:
	cd supplement; make

clean:
	rm *.o *-test
