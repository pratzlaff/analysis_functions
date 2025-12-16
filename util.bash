#! /bin/bash

i_of_n() {
    local i=$1
    local n=$2

    args=( "${@}" )
    unset args[0]
    unset args[1]
    args=( "${args[@]}" )
    argsstr=${args[@]}
    subset=$(perl -I/data/legs/rpete/flight/ARLac/src -MARLac=split_array -le "\$i=$i; \$n=$n; @args=qw!$argsstr!;"'print join " ", @{(split_array($n, @args))[$i-1]};')
    echo $subset
}
