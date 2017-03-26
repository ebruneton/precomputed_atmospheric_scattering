@echo off

REM create output dir
if not exist %~dp3 md %~dp3

REM remove comments
%1 -e "/^\/\*/,/\*\/$/d" -e "/^ *\/\//d" -e "/^$/d" %2 > %3.tmp0

REM for each line add quotes \r\n and \ : line -> "{line}\r\n" \
%1 "s/^\(.*\)$/\"\1\\r\\n\"\\/" %3.tmp0 > %3.tmp1

REM add variable declaration and ;
echo const char* %~n2_glsl = \> %3
type %3.tmp1 >> %3
echo ""; >> %3
