
generateProject and CMakelists.txt requires cmake, freeglut, glew and gnuwin32 in a folder in the current folder.
\precomputed_atmospheric_scattering\platform\windows\external\

ex:
external\cmake\bin\cmake.exe
external\freeglut\lib\x64\freeglut.lib
external\glew\lib\Release\x64\glew32.lib

for gnuwin32, only sed.exe is required:

external\gnuwin32\libiconv2.dll
external\gnuwin32\libintl3.dll
external\gnuwin32\regex2.dll
external\gnuwin32\sed.exe

tested with:

cmake 3.5
freeglut 3.0.0-1.mp
glew-2.0.0
sed-4.2.1
