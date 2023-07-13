#--------------------------------------- Select some defaults
CC       = gcc
CFLAGS   = -g -O0 -march=native -mavx -Wall -Wextra -std=c99

#--------------------------------------- Select target computer
SYSTYPE = "Accre"

#--------------------------------------- Adjust settings for target computer
ifeq ($(SYSTYPE),"Accre")
CC        = gcc
CFLAGS    = -O3 -march=native -mavx -Wall -Wextra -std=c99
LDFLAGS   = -lm -lqhull
QLIB      = -L$(QHULLSRCDIR)
QINC      = -I$(QHULLSRCDIR) 
endif


SRCS    = main.c io.c vozutil.c allvars.c
OBJS    = $(SRCS:.c=.o)
LIBOBJS = $(addprefix lib/,$(OBJS))
LIBS    = $(addprefix lib/lib,$(SRCS:.c=.so))
EXEC    = vorovol
INCL    = allvars.h proto.h Makefile $(QHULLSRCDIR)/qhull_a.h

#CFLAGS  = $(CFLAGS) #-fopenmp

#export LD_LIBRARY_PATH=$(QINC):$(LD_LIBRARY_PATH)

all: $(SRCS) $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(QLIB) $(QINC) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

%.o: %.c $(INCL)
	$(CC) $(QLIB) $(QINC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

lib/%.o: %.c $(INCL)
	$(CC) $(QLIB) $(QINC) -fpic $(CFLAGS) -c $< -o $@ $(LDFLAGS)

lib%.so: lib/%.o
	$(CC) $(QLIB) $(QINC) $(CFLAGS) $(LDFLAGS) -shared -o lib/$@ $< 

clean:
	rm -f $(OBJS) $(EXEC) $(LIBS) $(LIBOBJS)


# $(EXEC): $(OBJS)
#         $(CC) $(CFLAGS) $(QINC) $(OBJS) $(LDFLAGS)   -o  $(EXEC)  

#$(OBJS): $(INCL)



