#!/bin/sh

# NOTE: to use local API description, put gl.xml into glad/ directory so that glad sees it in its $PWD
cd "`dirname $0`/glad"
glad --out-path=. --generator=c --omit-khrplatform --api="gl=3.3" --profile=core --extensions=
mv src/glad.c src/glad.cc
