######################################################################################
# Makefile for building: CSPlibraries
# Usage:
# "make CPLEX" to produce CSPlibCPLEX
# "make XPRESS" to produce CSPlibXPRESS
# "make SCIP" to produce CSPlibSCIP
# "make" to produce all three libraries in one go
# use option "32BIT=true" to compile for 32 bit system
######################################################################################

####### Compiler, tools and options
# Environment
CC               = g++ -fPIC
CXX              = g++ -fPIC
WINDRES          = x86_64-w64-mingw32-windres
MKDIR            = mkdir
RM               = rm -f
CP               = cp -p

32BIT            = false
#32BIT            = false

ifeq ($(32BIT), false) # 64 bit assumed
    BITS         = -m64 -D_LP64
    ARCH         = x86_64
    CND_PLATFORM = MinGW-Windows64
    GNUDIR       = C:/Progra~1/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin
else
    BITS         = -m32
    ARCH         = x86
    CND_PLATFORM = MinGW-Windows
    GNUDIR       = C:/Progra~2/mingw-w64/i686-8.1.0-win32-sjlj-rt_v6-rev0/mingw32/bin
endif

# Macros
CND_DLIB_EXT     = so
CND_CONF         = Debug
CND_DISTDIR      = dist
CND_BUILDDIR     = build

# Object Directory
OBJECTDIR=$(CND_BUILDDIR)/$(CND_CONF)/$(CND_PLATFORM)
# for inclusion of ICallback.h and IProgressListener.h
INCTAUPATH       = -I../tauhitas/src

####### Object Files
OBJECTS          = $(OBJECTDIR)/src/Cspbridg.o $(OBJECTDIR)/src/Cspcard.o $(OBJECTDIR)/src/Cspmain.o \
                   $(OBJECTDIR)/src/Cspnet.o $(OBJECTDIR)/src/Jjsolver.o $(OBJECTDIR)/src/MT1RC.o \
                   $(OBJECTDIR)/src/My_time.o $(OBJECTDIR)/src/cspback.o $(OBJECTDIR)/src/cspbranc.o \
                   $(OBJECTDIR)/src/cspcapa.o $(OBJECTDIR)/src/cspcover.o $(OBJECTDIR)/src/cspdebug.o \
                   $(OBJECTDIR)/src/cspgomo.o $(OBJECTDIR)/src/cspheur.o $(OBJECTDIR)/src/cspprep.o \
                   $(OBJECTDIR)/src/cspprice.o $(OBJECTDIR)/src/cspsep.o $(OBJECTDIR)/src/cspsolve.o \
                   #$(OBJECTDIR)/src/Versioninfo.o

STAMP            =# -DSTAMP

CXXFLAGS	 = -g -O2 -Wall -DMICROSOFT2 $(BITS)
#CXXFLAGS         = -ggdb -O0 -Wall -DMICROSOFT2 $(BITS)

CSPCPX           = libCSP_CPLEX
CPXFLAGS         = $(CXXFLAGS) -DCPLEX7 -DBUILD_CPXSTATIC
ifeq ($(32BIT),false)
    CPXDIR       = ../Solvers/Cplex/Cplex125/Windows/64bits
    CPXLIBS      = -L$(CPXDIR) -lcplex125
else
    #CPXDIR       = ../Solvers/Cplex/Cplex75
    #CPXLIBS      = -L$(CPXDIR)/lib -lcplex75
    CPXDIR       = ../Solvers/Cplex/Cplex125/Windows/32bits
    CPXLIBS      = -L$(CPXDIR) -lcplex125
endif
CPXINC           = -I$(CPXDIR)/include/ilcplex

CSPXPR           = libCSP_XPRESS
XPRFLAGS         = $(CXXFLAGS) -DXPRESS_13
#XPRDIR           = ../Solvers/XPress/XPress_28/$(ARCH)
XPRDIR           = ../Solvers/XPress/XPress_19
XPRINC           = -I$(XPRDIR)
XPRLIBS          = -L$(XPRDIR) -lxprl -lxprs

CSPSCIP          = libCSP_SCIP
SCIPFLAGS        = $(CXXFLAGS) -DVSCIP -Dsoplex
DIRLPS           = ../Solvers/scip-3.1.1
DIRSOPLEX        = ../Solvers/soplex-2.0.1
SOPLEXLIB        = soplex-2.0.1.linux.$(ARCH).gnu.opt
NLPILIB          = nlpi.cppad-3.1.1.linux.$(ARCH).gnu.opt
SCIPLIB          = scip-3.1.1.linux.$(ARCH).gnu.opt
OBJSCIPLIB       = objscip-3.1.1.linux.$(ARCH).gnu.opt
LPISPXLIB        = lpispx-3.1.1.linux.$(ARCH).gnu.opt
SCIPINC          = -I$(DIRLPS)/src -I$(DIRSOPLEX)/src
SCIPLIBS         = -L$(DIRLPS)/lib -L$(DIRSOPLEX)/lib -L$(DIRLPS)/lib -l$(OBJSCIPLIB) -l$(SCIPLIB) -l$(NLPILIB) -l$(LPISPXLIB) -l$(SOPLEXLIB)

.PHONY: all clean

all : CPLEX XPRESS SCIP

CPLEX :
	$(MKDIR) -p $(OBJECTDIR)/src
	$(MKDIR) -p $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)
	#$(WINDRES) ./src/Versioninfo.rc $(CND_BUILDDIR)/$(CND_CONF)/$(CND_PLATFORM)/src/Versioninfo.o
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/Cspbridg.o src/Cspbridg.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/Cspcard.o src/Cspcard.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Cspmain.o src/Cspmain.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Cspnet.o src/Cspnet.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Jjsolver.o src/Jjsolver.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/MT1RC.o src/MT1RC.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/My_time.o src/My_time.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspback.o src/cspback.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/cspbranc.o src/cspbranc.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/cspcapa.o src/cspcapa.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/cspcover.o src/cspcover.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/cspdebug.o src/cspdebug.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspgomo.o src/cspgomo.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspheur.o src/cspheur.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/cspprep.o src/cspprep.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/cspprice.o src/cspprice.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) -o $(OBJECTDIR)/src/cspsep.o src/cspsep.c
	$(CXX) -c $(STAMP) $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspsolve.o src/cspsolve.c
	$(CXX) -o $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)/$(CSPCPX).$(CND_DLIB_EXT) $(OBJECTS) $(STAMP) $(CPXFLAGS) $(CPXLIBS) -shared #-static-libgcc -static-libstdc++
	$(CP) $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)/$(CSPCPX).$(CND_DLIB_EXT) ../tauargus

XPRESS :
	$(MKDIR) -p $(OBJECTDIR)/src
	$(MKDIR) -p $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)
	$(WINDRES) ./src/Versioninfo.rc $(CND_BUILDDIR)/$(CND_CONF)/$(CND_PLATFORM)/src/Versioninfo.o
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/Cspbridg.o src/Cspbridg.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/Cspcard.o src/Cspcard.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Cspmain.o src/Cspmain.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Cspnet.o src/Cspnet.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Jjsolver.o src/Jjsolver.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/MT1RC.o src/MT1RC.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/My_time.o src/My_time.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspback.o src/cspback.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/cspbranc.o src/cspbranc.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/cspcapa.o src/cspcapa.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/cspcover.o src/cspcover.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/cspdebug.o src/cspdebug.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspgomo.o src/cspgomo.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspheur.o src/cspheur.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/cspprep.o src/cspprep.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/cspprice.o src/cspprice.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o $(OBJECTDIR)/src/cspsep.o src/cspsep.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspsolve.o src/cspsolve.c
	$(CXX) -o $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)/$(CSPXPR).$(CND_DLIB_EXT) $(OBJECTS) $(XPRFLAGS) $(XPRLIBS) -shared #-static-libgcc -static-libstdc++
	$(CP) $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)/$(CSPXPR).$(CND_DLIB_EXT) ../tauargus

SCIP : 
	$(MKDIR) -p $(OBJECTDIR)/src
	$(MKDIR) -p $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)	
	$(WINDRES) ./src/Versioninfo.rc $(CND_BUILDDIR)/$(CND_CONF)/$(CND_PLATFORM)/src/Versioninfo.o
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/Cspbridg.o src/Cspbridg.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/Cspcard.o src/Cspcard.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Cspmain.o src/Cspmain.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Cspnet.o src/Cspnet.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/Jjsolver.o src/Jjsolver.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/MT1RC.o src/MT1RC.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/My_time.o src/My_time.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspback.o src/cspback.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/cspbranc.o src/cspbranc.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/cspcapa.o src/cspcapa.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/cspcover.o src/cspcover.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/cspdebug.o src/cspdebug.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspgomo.o src/cspgomo.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspheur.o src/cspheur.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/cspprep.o src/cspprep.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/cspprice.o src/cspprice.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o $(OBJECTDIR)/src/cspsep.o src/cspsep.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o $(OBJECTDIR)/src/cspsolve.o src/cspsolve.c
	$(CXX) -o $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)/$(CSPSCIP).$(CND_DLIB_EXT) $(OBJECTS) $(SCIPFLAGS) $(SCIPLIBS) -shared  -Wl,-Bsymbolic #-static-libgcc -static-libstdc++
	$(CP) $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)/$(CSPSCIP).$(CND_DLIB_EXT) ../tauargus
	
clean:
	$(RM) -r $(CND_BUILDDIR)/$(CND_CONF)
	$(RM) $(CND_DISTDIR)/$(CND_CONF)/$(CND_PLATFORM)/*.$(CND_DLIB_EXT)
