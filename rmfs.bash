#! /bin/bash

# for lsfparm_files
source /data/legs/rpete/flight/analysis_functions/caldb_files.bash

rmf_dir=/data/legs/rpete/flight/rmfs

generate_rmfs() {
    local pha2="$1"
    local lsfparm_dir="$2"
    local lsfparm_version="$3"

    local detnam=$(dmkeypar "$pha2" detnam echo+)

    local lsfparm_set=0
    case $# in
	1) ;;
	2|3) lsfparm_set=1 ;;
	*) echo 'Usage: generate_rmfs pha2 [lsfparm_directory [lsfparm_version]]' 1>&2
	   return 1
	   ;;
    esac

    local lsfparm_det
    local detsubsys
    local filename_det
    case "$detnam" in
	HRC-I) detsubsys=HRC-I ; lsfparm_det=hrc ; filename_det=HRC-I ;;
	HRC-S*) detsubsys=HRC-S2 ; lsfparm_det=hrc ; filename_det=HRC-S ;;
	ACIS-4*) detsubsys=ACIS-S3 ; lsfparm_det=acis ; filename_det=ACIS-S ;;
	*) echo "generate_rmfs() - unrecognized DETNAM='$detnam'" 1>&2
	   return 1 ;;
    esac

    local tg_m
    local tg_part
    local grating_arm
    local rmffile
    
    dmlist "$pha2"'[cols tg_m, tg_part]' data,raw | grep -v '^#' | while read line
    do
	read tg_m tg_part <<<$(echo "$line")
	case "$tg_part" in
	    1) grating_arm=HEG ;;
	    2) grating_arm=MEG ;;
	    3) grating_arm=LEG ;;
	    *) echo "generate_rmfs() - unrecognized TG_PART=$tg_part in $pha2" 1>&2
	       return 1 ;;
	esac

	#
	# Special cases
	#
	[ "$detnam" == HRC-I ] && {
	    [ "$grating_arm" != LEG ] && {
		detsubsys=ACIS-S3
		lsfparm_det=acis
	    }
	    [ "$grating_arm" == LEG ] && lsfparm_set=1
	}

	rmffile="${filename_det}-${grating_arm}_${tg_m}.rmf"

	#[ -f "$rmffile" ] && continue

	[ $lsfparm_set -ne 0 ] && {
	    local lsfparm_neg
	    local lsfparm_pos
	    read lsfparm_neg lsfparm_pos <<<$(lsfparm_files $lsfparm_det ${grating_arm} $lsfparm_dir $lsfparm_version)

	    local lsfparm="$lsfparm_pos"
	    [ $tg_m -lt 0 ] && lsfparm="$lsfparm_neg"

	    local lsf_string
	    case "$grating_arm" in
		HEG) lsf_string=HETG_0011 ;;
		MEG) lsf_string=HETG_1100 ;;
		LEG) lsf_string=LETG_1111 ;;
		*) return 1;; # should have errored out in earlier case "$tg_part"
	    esac

	    echo pset ardlib AXAF_${lsf_string}_LSF_FILE="$lsfparm"
	    pset ardlib AXAF_${lsf_string}_LSF_FILE="$lsfparm"
	}
	
	mkgrmf \
	    grating_arm="$grating_arm" \
	    order="$tg_m" \
            outfile="$rmffile" \
            srcid=1 \
            detsubsys=$detsubsys \
            threshold=1e-06 \
            obsfile="$pha2"'[SPECTRUM]' \
            regionfile="$pha2" \
            wvgrid_arf=compute \
            wvgrid_chan=compute \
            clobber=yes

    done
}
