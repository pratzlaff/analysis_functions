#! /bin/bash

tmppdir() {
    local dir='/tmp'
    [ -z "$1" ] || dir="$1"
    pdir=$(tmpdir "$dir")
    PFILES="$pdir;"$(echo "$PFILES" | perl -F';' -le 'print $F[1]')
    echo "$PFILES" 1>&2
    echo "$pdir"
}

# taken from https://stackoverflow.com/questions/4632028/how-to-create-a-temporary-directory

# deletes the temp directory
#cleanup() {      
#    rm -rf "$1"
#    echo "Deleted temp working directory $1" 1>&2
#}

tmpdir() {

    local DIR='/tmp'
    [ -z "$1" ] || DIR="$1"

    # the temp directory used, within $DIR omit the -p parameter to
    # create a temporal directory in the default location
    WORK_DIR=`mktemp -d -p "$DIR"`

    # check if tmp dir was created
    if [[ ! "$WORK_DIR" || ! -d "$WORK_DIR" ]]; then
	echo "Could not create temp dir"
	exit 1
    fi


    # register the cleanup function to be called on the EXIT signal
    #trap cleanup EXIT

    echo "$WORK_DIR"
}
