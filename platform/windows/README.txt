Cmake setup for windows

generate_project and CMakelists.txt requires cmake, freeglut, glew and gnuwin32 in precomputed_atmospheric_scattering\platform\windows\external\

ex:
precomputed_atmospheric_scattering\platform\windows\external\cmake\bin\cmake.exe
precomputed_atmospheric_scattering\platform\windows\external\freeglut\lib\x64\freeglut.lib
precomputed_atmospheric_scattering\platform\windows\external\glew\lib\Release\x64\glew32.lib

For gnuwin32, only sed.exe is required:

precomputed_atmospheric_scattering\platform\windows\external\gnuwin32\libiconv2.dll
precomputed_atmospheric_scattering\platform\windows\external\gnuwin32\libintl3.dll
precomputed_atmospheric_scattering\platform\windows\external\gnuwin32\regex2.dll
precomputed_atmospheric_scattering\platform\windows\external\gnuwin32\sed.exe

Tested with:

cmake 3.5
freeglut 3.0.0-1.mp
glew-2.0.0
sed-4.2.1
