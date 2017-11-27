#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo
    echo "        ./intall.sh -o destination_dir"
    echo
    echo "        Options :"
    echo "        -o: string   [the object directory of installing INTEGRATE-Neo]"
    echo
    exit
fi

args=("$@")

destination_dir=${args[1]}

echo "making destination dir..."

if [ ! -d $destination_dir/tmp ]; then
    mkdir -p $destination_dir"/tmp"
fi

echo "copying files..."

cp -r src/* $destination_dir
cd $destination_dir

echo "Compiling "

cd tmp
mv ../BedpeAnnotator ./
mkdir BedpeAnnotator_build
cd BedpeAnnotator_build
cmake ../BedpeAnnotator -DCMAKE_BUILD_TYPE=release
make

cd ..
cp BedpeAnnotator_build/bin/fusionBedpeAnnotator ../
cd ..

echo "adding +x mode"
chmod +x *py

echo "removing temporary files"

rm -rf tmp

echo 
echo
echo "Please make you you have installed the prerequisites. See details at https://github.com/ChrisMaherLab/INTEGRATE-Vis."
echo "Now enjoy!"
echo 
echo "You can run python "$destination_dir"/Integrate-vis.py" 
echo
