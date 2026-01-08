## Create WB polarimeric portraits with psradd (PSRCHIVE) or ppalin (PulsePortraiture; RECOMMENDED)
## Update "YOUR SET UP" below with your paths and datafiles and let it run...
##
## Run as:
## >> bash create_portrait.sh psr

psr=$1

#######################################
### YOUR SET UP 
#######################################

ext='ar'   # your data file extensions (e.g. obs1.ar)

datadir="your_data_dir"
resultdir="results/${psr}"

porttype="PP"  # "PSR" for PSRCHIVE or "PP" for ppalign
use_all=True   # if True, portrait will be generated from all available observations, if False will use $obs_num best observations (highest SNR)
obs_num=50

port_outfile="${resultdir}/${psr}.port_${porttype}"

#######################################
### PREPARE FOR ANALYSIS
#######################################

echo ""
echo "##########################################################"
echo "                 Procesing PSR ${psr}                    "
echo "##########################################################"
echo ""

mkdir -p ${resultdir}

metafile=${resultdir}/${psr}.meta #list of all observations
portmeta=${resultdir}/${psr}_port.meta #list of best $obs_num observations (ordered by SNR)
templ_ini=${resultdir}/${psr}.FT.sm #initial, smooth guess template (average profile)

touch ${metafile}
touch ${portmeta}

find "${datadir}" -type f -name "*.${ext}" ! -name "best_obs.${ext}" -print > "${metafile}"

#Create port metafile with obs_num
rm ${datadir}/best_obs.${ext}
psrstat -Q -j FTp -c snr ${datadir}/*.${ext} | sort -g -k2,2 |tail -n ${obs_num} |awk '{print $1}' > ${portmeta}

echo "Searching for highest SNR observation..."
# Choose 1 highest SNR observation
best=$(psrstat -Q -j DFTp -c snr ${datadir}/*.${ext} |sort -g -k2,2 |tail -n 1 | awk '{print $1}')
cp ${best} ${datadir}/best_obs.${ext}
echo ""

echo "Creating a smooth initial template..."
pam -e FT -FT ${best}
psrsmooth -W ${datadir}/*.FT
rm ${datadir}/*.FT
mv ${datadir}/*sm ${resultdir}/${psr}.FT.sm
echo ""

#######################################
### CREATE PORTRAIT
#######################################

if [[ "$use_all" == "true" ]]; then
    file_list=${metafile}
else
    file_list=${portmeta}
fi

if [[ "$porttype" == "PSR" ]]; then
    check_dmc=$(vap -c dmc "$(head -n 1 ${file_list})" | tail -n 1 | awk '{print $2}')
    if [ "$check_dmc" -eq 0 ]; then
        echo "Data is not dedispersed! I'll fix it."
	    pam -D -M ${file_list} -m
    fi

    psradd -TP -M ${file_list} -o ${port_outfile} -v

    #Convert to state="Stokes" in case data was in other state
    pam -S ${port_outfile} -m

elif [[ "$porttype" == "PP" ]]; then
    python3 ppalign_pol.py -M ${file_list} -I ${templ_ini} -C 15 -o ${port_outfile} --place 0.5 --poln --niter 3 --resultdir ${resultdir} 
    fi

echo ""
echo "##########################################################"
echo "#             Portrait created!                           "
echo "##########################################################"
echo ""

