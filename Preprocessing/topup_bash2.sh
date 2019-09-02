# !/bin/bash

export topupDIR=~/Desktop/TempMem_YJS/topup_data
#export allParticipants=('fmri101' 'fmri102' 'fmri103' 'fmri104' 'fmri105' 'fmri106' 'fmri107' 'fmri108' 'fmri109' 'fmri110' 'fmri111' 'fmri112' 'fmri113')

#export allParticipants=('fmri102' 'fmri103' 'fmri104' 'fmri105' 'fmri106' 'fmri107' 'fmri108' 'fmri109' 'fmri110' 'fmri111' 'fmri112' 'fmri113')
#export normalParticipants=('fmri102' 'fmri104' 'fmri105' 'fmri107' 'fmri108' 'fmri111' 'fmri112' 'fmri113')
export allParticipants=('fmri101' 'fmri103' 'fmri106' 'fmri109' 'fmri110')
export runs=('PRE_REST' 'ENC-01' 'ENC-02' 'ENC-03' 'ENC-04' 'ENC-05' 'ENC-06' 'POST_REST')
export ind=('topup1' 'topup2')
for isub in {0..4}
do
    ct=0
    for j in {0..1}
    do
    # go to each subject's folder
    export curr_dir=$topupDIR/${allParticipants[isub]}/day1/encoding/${ind[j]}
    cd $curr_dir


    # if-statement to differentiate subjects with diff. sequences
#if [[ ${doubleParticipants[${ct}]} == ${allParticipants[isub]} ]]; then



    # if participants are double
#   else
    echo "start topup correction"
    echo "now working on " ${allParticipants[isub]}



    fslmerge -t se_epi_merged *FM_P* *FM_A* #fmap_PA fmap_AP
    topup --imain=$curr_dir/se_epi_merged.nii.gz --datain=$topupDIR/parameter.txt --config=b02b0.cnf --out=topup_result --iout=topup_result_example
    echo "topup correction was done"
    if [ ${j} == 0 ]; then
        if [ ${isub} == 0 ]; then
            stoppoint=1
        elif [ ${isub} == 1 ]; then
            stoppoint=0
        elif [ ${isub} == 2 ]; then
            stoppoint=1
        elif [ ${isub} == 3 ]; then
            stoppoint=1
        elif [ ${isub} == 4 ]; then
            stoppoint=3
        fi
    elif [ ${j} == 1 ]; then
        if [ ${isub} == 0 ]; then
            stoppoint=5
        elif [ ${isub} == 1 ]; then
            stoppoint=6
        elif [ ${isub} == 2 ]; then
            stoppoint=5
        elif [ ${isub} == 3 ]; then
            stoppoint=5
        elif [ ${isub} == 4 ]; then
            stoppoint=3
        fi

    fi



        for irun in $(seq 0 $stoppoint)#{0..$stoppoint}
        do
        # if participants are normal
        echo "now working on " ${allParticipants[isub]} " and " ${runs[ct]}
        applytopup --imain=`ls $curr_dir/*${runs[ct]}*` --datain=$topupDIR/parameter.txt --inindex=2 --topup=topup_result --method=jac --out=${runs[ct]}
        echo "move to next RUN " ${runs[$((ct+1))]}
        ct=$((ct+1))
        done
#   fi
    done
done
