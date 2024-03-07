@echo off

setlocal
if not exist "build_cmake" mkdir "build_cmake"
cmake ^
   -D CMAKE_BUILD_TYPE=Debug ^
   -D BUILD_SHARED_LIBS=OFF ^
   -D enable-openmp=OFF ^
   -D enable-sym-ev-routine=OFF ^
   -D enable-single=OFF ^
   -D enable-gpu-code=OFF ^
   -D enable-cuda=ON ^
   -D enable-gpu-double=OFF ^
   -S . -B build_cmake %~1

set "cleanfirst=--clean-first"
REM set "cleanfirst="

cmake --build build_cmake --config Debug %cleanfirst%
cmake --build build_cmake --config Release %cleanfirst%
endlocal
