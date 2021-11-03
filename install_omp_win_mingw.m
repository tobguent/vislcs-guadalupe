% The scripts use MEX, which has to be installed once with the command
%   mex -setup C++
% Afterwards, this file can be executed to compile the code.
% -----------------------------------------------------------------------------------
% Use this script if you use a 'MinGW64' compiler on Windows and want OpenMP support.
% -----------------------------------------------------------------------------------
% Note that you will have to adjust the path to your 'libgomp.a' files!
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cstreamline.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cevstreamline.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cpathline.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cstreakline.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/ctimeline.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/clavd.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/civd.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cflowmap.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cftle.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/cvorticity.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/clic.cpp;
mex C:\ProgramData\MATLAB\SupportPackages\R2021a\3P.instrset\mingw_w64.instrset\lib\gcc\x86_64-w64-mingw32\6.3.0\libgomp.a CFLAGS='$CFLAGS -fopenmp' -LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" ./cpp/ctexadvect.cpp;

disp('compiling finished');