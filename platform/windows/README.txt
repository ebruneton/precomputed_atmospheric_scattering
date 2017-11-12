Cmake setup for Visual Studio.
Tested on Windows 10 with VS2017 x64.

The required dependencies are cmake, freeglut and sed.
A script is provided to download the versions of those dependencies available at the time when this is published.

They need to be placed in "precomputed_atmospheric_scattering\platform\windows\external"
following this hierarchy:

precomputed_atmospheric_scattering\platform\windows\external\cmake\bin\cmake.exe
precomputed_atmospheric_scattering\platform\windows\external\freeglut\lib\x64\freeglut.lib
precomputed_atmospheric_scattering\platform\windows\external\sed\bin\sed.exe

Tested with:

cmake 3.9
freeglut 3.0.0-1.mp
sed-4.2.1

download_build_run.bat: download dependencies, generate project, build and run.

download_dependencies.bat: download the dependencies. 
it's using some Powershell modules which might only be available on Windows 8 and above.

generate_project.bat: generate Visual Studio project, requires the dependencies.

build.bat: build the project in Release using the cmake project, requires generate_project.

run.bat: run the built project.


