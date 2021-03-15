
cd ${1}
for i in $(seq 4); do
	cd $i
	touch name.txt
	cd ..
done
cd ..
