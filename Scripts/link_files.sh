#!/bin/bash
#May 30, 2017 by Teresita M. Porter
#script to create a directory of files using symbolic links only

echo "Please enter path containing original files (OMIT final '/'):"
read path

cwd=$(pwd)
echo $cwd

for f in $path/*.gz
	do
		base=${f//$path\//}
		ln -s $path/$base $cwd/$base
	done
