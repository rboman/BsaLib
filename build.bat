@echo off

setlocal
if not exist "build_cmake" mkdir "build_cmake"
cmake ^
   -D enable-openmp=ON ^
   -D enable-sym-ev-routine=OFF ^
   -D enable-single=OFF ^
   -D enable-gpu-code=OFF ^
   -D enable-cuda=ON ^
   -D enable-gpu-double=OFF ^
   -S . -B build_cmake %~1
endlocal
