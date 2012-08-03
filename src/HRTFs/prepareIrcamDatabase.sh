top=1018
mkdir -p IRCAM
(
cd IRCAM
for ((i=1016;i<$top;i++)); do
	wget -c 'ftp://ftp.ircam.fr/pub/IRCAM/equipes/salles/listen/archive/SUBJECTS/IRC_'$i'.zip';
	if [ ! -e IRC_$i''.zip ]
	then
		echo -e "\033[31mError downloading IRC_$i.zip \033[0m" >&2
		continue
	fi
	unzip IRC_$i''.zip
	mkdir -p IRCAM_$i''_{L,R}
	false > ../ircam${i}L.hrtfs
	false > ../ircam${i}R.hrtfs
	for a in COMPENSATED/WAV/IRC_$i''_C/IRC_$i''_C_R0195_T*
	do
		base=`basename $a`
		elevat=`echo $base | sed "s%.*_P\\([0-9]*\\).wav%\1%"`
		azimuth=`echo $base | sed "s%.*_T\\([0-9]*\\).*%\1%"`
		sox $a -c 1 -2 IRCAM_$i''_L/$base mixer -l
		sox $a -c 1 -2 IRCAM_$i''_R/$base mixer -r
		echo $elevat $azimuth IRCAM/IRCAM_$i''_L/$base >> ../ircam${i}L.hrtfs 
		azimuth_without_zero=`echo $azimuth | sed 's/^0//'`
		echo $elevat $((360-azimuth_without_zero)) IRCAM/IRCAM_$i''_R/$base >> ../ircam${i}R.hrtfs 
	done
	rm -rf RAW COMPENSATED
done
)


