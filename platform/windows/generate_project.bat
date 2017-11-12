@echo off

set CURDIR=%~dp0
set CMAKEPATH=%CURDIR%external\cmake\bin
set CMAKEFOLDER=_intermediate

if exist %CMAKEFOLDER% rmdir %CMAKEFOLDER% /s /q
mkdir %CMAKEFOLDER%

pushd %CMAKEFOLDER%

%CMAKEPATH%\cmake -G"Visual Studio 15 2017 Win64" ..

popd

if NOT '%1' == 'NOPAUSE' pause