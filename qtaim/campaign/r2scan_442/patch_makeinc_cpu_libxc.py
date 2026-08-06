#!/usr/bin/env python3
"""Transform the QE GPU make.inc into a CPU + libxc make.inc.
Strips CUDA (__CUDA, CUDA_* flags/libs, -cuda/-gpu/-acc), adds -D__LIBXC +
libxc include/libs (libxcf03 = the F2003 module 'xc_f03_lib_m' QE XClib USEs),
points BLAS/LAPACK at NVHPC CPU compilers/lib. Keeps nvfortran/mpif90 toolchain."""
import sys, re
src = sys.argv[1]   # /workspace/qe-7.5-libxc/make.inc
dst = sys.argv[2]   # overwrite same
s = open(src).read()
NVLIB = "/workspace/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/lib"

def sub(line, old, new):
    assert old in line, f"NOT FOUND in line: {old!r} | line={line!r}"
    return line.replace(old, new)

out = []
for line in s.splitlines(keepends=True):
    L = line
    if L.startswith("DFLAGS"):
        L = sub(L, "-D__PGI -D__CUDA -D__FFTW -D__MPI -D__MPI_MODULE",
                   "-D__PGI -D__FFTW -D__MPI -D__MPI_MODULE -D__LIBXC")
    elif L.startswith("IFLAGS"):
        L = sub(L, "-I. -I$(TOPDIR)/include",
                   "-I. -I$(TOPDIR)/include -I/workspace/libxc/include")
    elif L.startswith("CUDA_F90FLAGS"):
        L = "CUDA_F90FLAGS =\n"
    elif L.startswith("CUDA_CFLAGS"):
        L = "CUDA_CFLAGS =\n"
    elif L.startswith("CFLAGS"):
        L = sub(L, " $(CUDA_CFLAGS)", "")
    elif L.startswith("F90FLAGS"):
        L = sub(L, " $(CUDA_F90FLAGS)", "")
    elif L.startswith("LDFLAGS"):
        L = "LDFLAGS        =  -mp\n"
    elif L.startswith("BLAS_LIBS"):
        L = f"BLAS_LIBS      =  -L{NVLIB} -lblas \n"
    elif L.startswith("LAPACK_LIBS"):
        L = f"LAPACK_LIBS = -L{NVLIB} -llapack -lblas \n"
    elif L.startswith("CUDA_LIBS"):
        L = "CUDA_LIBS =\n"
    elif L.startswith("CUDA_EXTLIBS"):
        L = "CUDA_EXTLIBS =\n"
    elif L.startswith("LIBXC_LIBS"):
        L = "LIBXC_LIBS     = -L/workspace/libxc/lib -lxcf03 -lxc\n"
    out.append(L)
open(dst, "w").write("".join(out))
print("patched make.inc -> CPU + libxc")
# echo the critical lines for verification
for key in ("DFLAGS","IFLAGS","F90FLAGS","CFLAGS","LDFLAGS","LIBXC_LIBS","BLAS_LIBS","LAPACK_LIBS","CUDA_LIBS","CUDA_F90FLAGS"):
    for l in "".join(out).splitlines():
        if l.startswith(key+" ") or l.startswith(key+"="):
            print("  ", l); break
