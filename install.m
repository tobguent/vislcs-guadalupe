% The scripts use MEX, which has to be installed once with the command
%   mex -setup C++
% Afterwards, this file can be executed to compile the code.
% ---------------------------------------------------------------------
% This is the basic compiler script, which does not use OpenMP support.
% ---------------------------------------------------------------------
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/cstreamline.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/cevstreamline.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/cpathline.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/cstreakline.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/ctimeline.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/clavd.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/civd.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/cflowmap.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/cftle.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/cvorticity.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/clic.cpp;
mex COPTIMFLAGS="-O3 -DNDEBUG" ./cpp/ctexadvect.cpp;

disp('compiling finished');