@echo off

setlocal
if not exist "build_cmake" mkdir "build_cmake"
cmake ^
   -D enable-sym-ev-routine=OFF ^
   -D enable-gpu-code=ON ^
   -S . -B build_cmake %~1
endlocal
