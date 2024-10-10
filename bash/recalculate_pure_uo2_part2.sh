cd /media/jarinf/Research1/uo2/grain_growth/cylindrical
for i in $(fd find_grains_input2.txt); do 
    (
    cd $(dirname ${i})
    n_fields=$(head -n1 find_grains_input2.txt | awk '{print NF}')
    if (( n_fields == 6 )); then
        awk 'NR==1 {for (i=2;i<NF;i++) printf $i " "; print $NF} NR > 1 {print $0}' find_grains_input2.txt > tmp
        mv tmp find_grains_input2.txt
        find_grains_v2 find_grains_input2.txt -i 2 -e 0 -n 0.0
        pot=$(echo ${i} | awk -F'/' '{print $3}')
        if [[ ${pot} == "Basak" ]]; then
            pnum=1
        elif [[ ${pot} == "Cooper" ]]; then
            pnum=17
        else
            echo "Unknown potential ${pot} in ${i}"
            continue
        fi
        calculate_grain_area grain_area_input.txt -p ${pnum}
        mobility_plots.py area_data.txt --simple -F
    fi
    )
done

cd ${OLDPWD}