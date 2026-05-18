#! /bin/bash

# for $rmf_dir
source /data/legs/rpete/flight/analysis_functions/rmfs.bash

asphist_dir=./asphists

generate_garfs() {
    local indir="$1"
    local outdir="$2"

    mkdir -p "$outdir"

    local evt2=$(\ls "$indir"/*_evt2.fits)
    local pha2=$(\ls "$indir"/*_pha2.fits)
    local asol=$(asol_stack "$indir"/../primary)

    local dtf1

    local obsid=$(dmkeypar "$pha2" obs_id echo+)
    local grating=$(dmkeypar "$pha2" grating echo+)

    local detnam=$(dmkeypar "$pha2" detnam echo+)
    case "$detnam" in
	HRC-*) dtf1=$(\ls "$indir"/../primary/*_dtf1.fits* | tail -1) ;;
	ACIS-[456789]*) detnam=ACIS-S ; dtf1="$evt2" ;;
	*) echo "generate_garfs: unrecognized DETNAM=$detnam" 1>&2 ; return 1 ;;
    esac

    local row=0
    dmlist "$pha2"'[cols tg_m, tg_part]' data,raw | \grep -v '^#' | while read line
    do
	(( row+=1 ))
	read tg_m tg_part <<<$(echo "$line")

	[ ! -z "$minorder" ] && [ ${tg_m#-} -lt $minorder ] && continue
	[ ! -z "$maxorder" ] && [ ${tg_m#-} -gt $maxorder ] && continue

	case "$tg_part" in
	    1) grating_arm=HEG ;;
	    2) grating_arm=MEG ;;
	    3) grating_arm=LEG ;;
	    *) echo "generate_garfs: unrecognized TG_PART=$tg_part in $pha2" 1>&2 ; return 1 ;;
	esac

	local rmf="${rmf_dir}/${detnam}-${grating_arm}_${tg_m}.rmf"

	local bpix1
	case "${detnam,,}" in
	    hrc*)  bpix1=$(\ls "$indir"/*_bpix1.fits 2>/dev/null)
		   [ -z "$bpix1" ] && bpix1=$(\ls "$indir"/../secondary/*_bpix1.fits* | tail -1)
		   echo pset ardlib AXAF_${detnam}_BADPIX_FILE="${bpix1}[BADPIX]"
		   pset ardlib AXAF_${detnam}_BADPIX_FILE="${bpix1}[BADPIX]"
		   ;;
	    acis*) bpix1="$indir/acis_repro_bpix1.fits"
		   ;;
	esac

	echo bpix1=$bpix1

	punlearn mkgarf
	punlearn fullgarf
	[ ! -z "$ardlibqual" ] && pset fullgarf ardlibqual="$ardlibqual"
	rm -f "$outdir/$obsid"_${grating_arm}_${tg_m}_garf.fits
	rm -f "$outdir/$obsid"__S?_${grating_arm}_${tg_m}.fits
	fullgarf \
            "$pha2" \
            "$row" \
            "$evt2" \
            "$asol" \
            "grid($rmf[cols ENERG_LO,ENERG_HI])" \
            "$dtf1" \
            "$bpix1" \
            "$outdir/${obsid}_" \
            maskfile=NONE \
            clobber=no

	# currently there's a CIAO bug which causes dmaddarf to not
	# retain the tg_m keyword when doing a copy of a single file
	# (e.g., HRC-I fullgarf output)
	# Update: as of CIAO 4.14 this appears to have been fixed
	false && dmhedit \
	    "$outdir/${obsid}_${grating_arm}_${tg_m}_garf.fits" \
	    filelist=none \
	    operation=add \
	    key=TG_M \
	    value=$tg_m \
	    datatype=indef

	rm -f "$outdir/$obsid"__S?_${grating_arm}_${tg_m}.fits
    done
}

generate_0th_arf() {
    local indir="$1"
    local outdir="$2"

    mkdir -p "$asphist_dir" "$outdir"

    local evt2=$(\ls "$indir"/*_evt2.fits)
    local pha2=$(\ls "$indir"/*_pha2.fits)

    local obsid=$(dmkeypar "$pha2" obs_id echo+)
    local detnam=$(dmkeypar "$pha2" detnam echo+)
    local grating=$(dmkeypar "$pha2" grating echo+)

    local detsubsys="$detnam"
    case "$detnam" in
	HRC-I) ;;
	HRC-S) detsubsys=HRC-S2 ;;
	ACIS-4*) detsubsys=ACIS-S3 ;;
	*) echo "generate_0th_arf: unrecognized DETNAM=$detnam" 1>&2
	   return 1;
	   ;;
    esac

    local x y
    read x y <<<$( dmlist "$evt2"'[region][col sky]' data,raw | head -2 | tail -1 )

    local rmffile="$rmf_dir/${detnam}-LEG_1.rmf"
    local arf="$outdir/${obsid}_0th.arf"
    local asphist="${asphist_dir}/${obsid}_asphist.fits"
    local dtf=$(\ls "$indir"/../primary/*_dtf1.fits*)

    punlearn asphist
    asphist \
        $(asol_stack "$indir"/../primary) \
        "$asphist" \
        "$evt2" \
        "$dtf"

    punlearn mkarf
    mkarf \
        asphistfile="$asphist" \
        outfile="$arf" \
        sourcepixelx="$x" \
        sourcepixely="$y" \
        engrid='grid('"$rmffile"'[cols ENERG_LO,ENERG_HI])' \
        obsfile="$pha2"'[SPECTRUM]' \
        detsubsys="$detsubsys" \
        grating=$grating \
        maskfile=NONE \
        verbose=0 \
        clobber+
}
