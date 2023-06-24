#! /bin/bash

centroid_sky() {
  infile="$1"
  x="$2"
  y="$3"
  r="$4"
  for i in {1..10}
    do
    read -r x y <<<$(dmstat "$infile""[(x,y)=circle($x, $y ,$r)][col x,y]" | grep mean | perl -anle 'print $F[2], "\t", $F[3]')
    echo $x $y
  done
}

centroid_eqpos() {
  infile="$1"
  x="$2"
  y="$3"
  r="$4"
  for i in {1..10}
    do
    read -r x y <<<$(dmstat "$infile""[(ra,dec)=circle($x, $y ,$r)][col ra, dec]" | grep mean | perl -anle 'print $F[2], "\t", $F[3]')
    echo $x $y
  done
}

