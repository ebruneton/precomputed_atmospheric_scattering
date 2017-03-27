@echo off

set SED=%~dp0external\sed\bin\sed.exe

REM create output dir
if not exist %~dp2 md %~dp2

REM remove comments
%SED% -e "/^\/\*/,/\*\/$/d" -e "/^ *\/\//d" -e "/^$/d" %1 > %2.tmp0

REM for each line add quotes \r\n and \ : line -> "{line}\r\n" \
%SED% "s/^\(.*\)$/\"\1\\r\\n\"\\/" %2.tmp0 > %2.tmp1

REM add variable declaration and ;
echo const char* %~n1_glsl = \> %2
type %2.tmp1 >> %2
echo ""; >> %2
