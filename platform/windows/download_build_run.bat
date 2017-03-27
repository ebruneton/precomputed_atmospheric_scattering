@echo off

call "download_dependencies.bat" NOPAUSE
call "generate_project.bat" NOPAUSE
call "build.bat" NOPAUSE
call "run.bat" NOPAUSE