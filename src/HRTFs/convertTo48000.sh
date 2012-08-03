echo ========== Duplicating the HRTF database =================
mkdir -p v48000
cp -r kreuzer2 v48000/.
cp -r IRCAM v48000/.
cp -r MIT_KEMAR v48000/.
cp *.hrtfs v48000/.

echo ========== Resampling to 48K =================
(cd v48000 && find -iname "*.wav" | while read a; do sox $a $a.temp.wav rate -v 48000; mv $a.temp.wav $a; done)

