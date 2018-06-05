#!/bin/bash
head -2 $1 | tail -1 > oneline
filesize=$(du -b $1 | cut -f -1)
linesize=$(du -b oneline | cut -f -1)
rm oneline
echo $(expr $filesize / $linesize)
