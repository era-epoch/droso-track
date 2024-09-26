#!/bin/bash
javac wrMTrck_.java -d build
cp -r build /mnt/c/users/jacob/documents/Flies/wrmTrackMod/
rm -r /mnt/c/users/jacob/documents/Flies/wrmTrackMod/eraMod
mv /mnt/c/users/jacob/documents/Flies/wrmTrackMod/build /mnt/c/users/jacob/documents/Flies/wrmTrackMod/eraMod