% The scripts use MEX, which has to be installed once with the command
%  mex -setup C++

% Afterwards, this file can be executed to compile the code.
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cstreamline.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cpathline.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cstreakline.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/ctimeline.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/clavd.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/civd.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cftle.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cvorticity.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/clic.cpp;
mex COMPFLAGS="$COMPFLAGS /openmp" CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/ctexadvect.cpp;
disp('compiling finished');
