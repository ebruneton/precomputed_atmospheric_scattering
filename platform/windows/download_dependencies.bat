@echo off
REM     https://cmake.org/files/v3.9/cmake-3.9.6-win64-x64.zip -> cmake
REM     http://sourceforge.net/projects/gnuwin32/files//sed/4.2.1/sed-4.2.1-bin.zip/download -> sed
REM     http://sourceforge.net/projects/gnuwin32/files//sed/4.2.1/sed-4.2.1-dep.zip/download -> sed
REM     http://files.transmissionzero.co.uk/software/development/GLUT/freeglut-MSVC.zip -> freeglut

set EXT=%~dp0external

echo Installing dependencies into %EXT%
if not exist %EXT% md %EXT%

REM ------------------ freeglut

set SRC=http://files.transmissionzero.co.uk/software/development/GLUT/freeglut-MSVC.zip
set ZIP=%EXT%\freeglut-MSVC.zip
set DST=%EXT%\.

if not exist %EXT%\freeglut (
echo Downloading freeglut from %SRC%...
powershell -Command "Start-BitsTransfer '%SRC%' '%ZIP%'"
powershell -Command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('%ZIP%', '%DST%'); }"
del %ZIP%
) else echo freeglut detected. skipping.

REM ------------------ sed

set SRC=http://sourceforge.net/projects/gnuwin32/files//sed/4.2.1/sed-4.2.1-bin.zip/download
set ZIP=%EXT%\sed-4.2.1-bin.zip
set DST=%EXT%\sed

if not exist %EXT%\sed\bin\sed.exe (
echo Downloading sed from %SRC%...
powershell -Command "Start-BitsTransfer '%SRC%' '%ZIP%'"
powershell -Command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('%ZIP%', '%DST%'); }"
del %ZIP%
) else echo sed.exe detected. skipping.

REM ------------------ sed-dep

set SRC=http://sourceforge.net/projects/gnuwin32/files//sed/4.2.1/sed-4.2.1-dep.zip/download
set ZIP=%EXT%\sed-4.2.1-dep.zip
set DST=%EXT%\sed

if not exist %EXT%\sed\bin\libiconv2.dll (
echo Downloading sed-dep from %SRC%...
powershell -Command "Start-BitsTransfer '%SRC%' '%ZIP%'"
powershell -Command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('%ZIP%', '%DST%'); }"
del %ZIP%
) else echo sed-dep detected. skipping.

REM ------------------ cmake

set SRC=https://cmake.org/files/v3.9/cmake-3.9.6-win64-x64.zip
set ZIP=%EXT%\cmake-3.9.6-win64-x64.zip
set DST=%EXT%\.

if not exist %EXT%\cmake (
echo Downloading cmake from %SRC%...
powershell -Command "Start-BitsTransfer '%SRC%' '%ZIP%'"
powershell -Command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('%ZIP%', '%DST%'); }"
rename %EXT%\cmake-3.9.6-win64-x64 cmake
del %ZIP%
) else echo cmake detected. skipping.

:end
if NOT '%1' == 'NOPAUSE' pause

