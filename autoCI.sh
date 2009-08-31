#!/bin/bash
# used to expedite RCS check in
# dearl@soe 19 May 2009
comment=$1

if [ ! "$comment" ]; then
    echo -e "Usage: autoCI.sh 'I changed all of foo to bar!'"
    exit 1
fi

ci -l -m"$comment" autoCI.sh
ci -l -m"$comment" fileHandling.c
ci -l -m"$comment" hashListFunctions.c
ci -l -m"$comment" linearAlg.c
ci -l -m"$comment" Makefile
ci -l -m"$comment" spia.c
ci -l -m"$comment" spia.h
ci -l -m"$comment" probFunctions.c

exit 0
