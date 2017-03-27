@echo off

set CURDIR=%~dp0
set CMAKEPATH=%CURDIR%external\cmake\bin
set CMAKEFOLDER=_intermediate

pushd %CMAKEFOLDER%

%CMAKEPATH%\cmake --build . --config Release

popd

if NOT '%1' == 'NOPAUSE' pause