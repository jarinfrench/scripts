#!/usr/bin/env bash

shopt -s extglob
LAMMPS_INPUT_FILE=print_potential_function_tmp.in

pair_style_sub_choice() {
  if [ -z "$1" ]; then
    >&2 echo "No parent potential specified. Exiting..."
    exit 1
  fi

declare -a sub_style_choices
  case $1 in
    adp) sub_style_choices+=("adp")
      ;;
    born) sub_style_choices+=("born" "born/coul/long" "born/coul/msm")
      ;;
    buck) sub_style_choices+=("buck" "buck/coul/cut" "buck/coul/long" "buck/coul/msm")
      ;;
    comb) sub_style_choices+=("comb" "comb3")
      ;;
    coul) sub_style_choices+=("coul/cut" "coul/long")
      ;;
    eam) sub_style_choices+=("eam" "eam/alloy" "eam/cd" "eam/fs")
      ;;
    lj) sub_style_choices+=("lj/cut" "lj/cut/coul/cut" "lj/cut/coul/long")
      ;;
    meam) sub_style_choices+=("meam/c")
      ;;
    morse) sub_style_choices+=("morse")
      ;;
    *) >&2 echo "Parent potential style $1 not implemented. Exiting..."
      exit 2
      ;;
  esac

  if [[ "${#sub_style_choices[@]}" -eq 1 ]]; then
    >&1 echo "${sub_style_choices[0]}"
  else
    PS3="Select the $1 sub style to use: "
    select sub_style in "${sub_style_choices[@]}"; do break; done
    >&1 echo ${sub_style}
  fi
}

echo "Select the units:"
select units in lj real metal si cgs electron micro nano; do break; done
echo "units ${units}" > "${LAMMPS_INPUT_FILE}"

echo "Select the atom style:"
select atom_style in angle atomic body bond charge dipole dpd \
  edpd mdpd tdpd electron ellipsoid full line meso molecular peri smd \
  sphere spin tri template hybrid; do break; done
echo "atom_style ${atom_style}" >> "${LAMMPS_INPUT_FILE}"
echo "dimension 3" >> "${LAMMPS_INPUT_FILE}"
echo "boundary p p p" >> "${LAMMPS_INPUT_FILE}"

echo "Select the (base) crystal structure:"
select crystal_structure in sc bcc fcc hcp diamond; do break; done

read -rp "What is the lattice parameter? " a0
while [[ ! ${a0} =~ ^[+]?[0-9]+\.?[0-9]*$ ]]; do
  read -rp "Please enter a valid positive number: " a0
done
echo "lattice ${crystal_structure} ${a0}" >> "${LAMMPS_INPUT_FILE}"
echo "region 1 block 0.0 10.0 0.0 10.0 0.0 10.0" >> "${LAMMPS_INPUT_FILE}"

read -rp  "How many atom types are in the simulation? " n_types
while [[ ! ${n_types} =~ ^[+]?[0-9]+$ ]]; do
  read -rp "Please enter a valid positive integer: " n_types
done
echo "create_box ${n_types} 1" >> "${LAMMPS_INPUT_FILE}"
echo "create_atoms 1 region 1" >> "${LAMMPS_INPUT_FILE}"

if [[ "${n_types}" -gt 1 ]]; then
  for i in $(seq 2 ${n_types}); do
    # lattice commands
    dialog --inputbox "Enter the lattice commands for atom type ${i}" 10 80 2>>${LAMMPS_INPUT_FILE}
    # create_atoms commands
    echo "create_atoms ${i} region 1" > ${LAMMPS_INPUT_FILE}
  done
fi

if [[ ${atom_style} == "charge" ]]; then
  while {
    read -a charges -rp "Enter the charges of the ${n_types} atom types (in order), space separated: "
    # Check if we have the right number of values
    while [[ "${#charges[@]}" -ne "${n_types}" ]]; do
      read -a charges -rp "Please enter ${n_types} numbers, separated by a space: "
    done

    # Check that each value is a valid number
    for i in "${charges[@]}"; do
      if [[ ! "${i}" =~ ^[+-]?[0-9]+\.?[0-9]*$ ]]; then
        echo -e "Please enter ${n_types} \e[3mnumbers\e[0m"
        break 1
      else
        break 2
      fi
    done
  } do true
  done

fi


echo "Select the pair potential style:"
potential_styles="adp born buck comb coul eam lj meam morse"
potential_styles_case=$(echo "@(${potential_styles// /|})")
potential_styles_array=("adp" "born" "buck" "comb" "coul" "eam" "lj" "meam" "morse")
pair_style=""

select pair_style_choice in hybrid ${potential_styles}; do
  #shellcheck disable=SC2254
  case ${pair_style_choice} in
    hybrid)
      pair_style="hybrid/overlay"
      cmd=(dialog --separate-output --checklist "Select styles:" 22 76 12)
      options=(1 "adp"   off
               2 "born"  off
               3 "buck"  off
               4 "comb"  off
               5 "coul"  off
               6 "eam"   off
               7 "lj"    off
               8 "meam"  off
               9 "morse" off)
      choices=$("${cmd[@]}" "${options[@]}" 2>&1 >/dev/tty)
      clear
      for choice in ${choices}; do
        pair_style="${pair_style} $(pair_style_sub_choice ${potential_styles_array[$((choice-1))]})"
      done
      ;;
    ${potential_styles_case})
      pair_style=$(pair_style_sub_choice ${pair_style_choice})
      ;;
  esac
  break
done

pair_styles_with_no_cutoff="@(adp|comb|eam|eam_alloy|eam_cd|eam_fs|meam_c)"
pair_styles_with_one_cutoff="@(born|buck|coul_cut|coul_long|lj_cut|morse)"
pair_styles_with_one_cutoff_and_one_optional_cutoff="@(born_coul_long|born_coul_msm|buck_coul_cut|buck_coul_long|buck_coul_msm|lj_cut_coul_cut|lj_cut_coul_long)"

cutoffs1="" # for the required cutoff
cutoffs2="" # for the optional cutoff
for style in ${pair_style}; do
  # ask for additional parameters in the pair_style command
  #shellcheck disable=SC2254
  case ${style//'/'/'_'} in
    hybrid_overlay) continue
      ;;
    comb3)
      read -rp "Include atomic polarization (Y|n)? " comb3_polar
      while [[ ! ${comb3_polar} =~ ^[YyNn]$ ]]; do
        read -rp "Unrecognized response. Include atomic polarization (Y|n)? " comb3_polar
      done
      ;;
    ${pair_styles_with_no_cutoff})
      cutoffs1="${cutoffs1} na"
      cutoffs2="${cutoffs2} na"
      ;;
    ${pair_styles_with_one_cutoff})
      read -rp "Enter the ${style} cutoff in ${units} units: " tmp_cutoff
      while [[ ! ${tmp_cutoff} =~ ^[+]?[0-9]+\.?[0-9]*$ ]]; do
        read -rp "Please enter a valid positive value: " tmp_cutoff
      done
      cutoffs1="${cutoffs1} ${tmp_cutoff}"
      cutoffs2="${cutoffs2} na"
      ;;
    ${pair_styles_with_one_cutoff_and_one_optional_cutoff})
      read -rp "Enter ${style} cutoff1 (and optional cutoff2) in ${units} units: " reqd_cut coul_cut
      while [[ ! ${reqd_cut} =~ ^[+]?[0-9]+\.?[0-9]*$ ]] || ([[ -n ${coul_cut} ]] && [[ ! ${coul_cut} =~ ^[+]?[0-9]+\.?[0-9]*$ ]]);
      do
        read -rp "Please enter valid positive value(s): " reqd_cut coul_cut
      done
      cutoffs1="${cutoffs1} ${reqd_cut}"
      if [[ -z ${coul_cut} ]]; then
        cutoffs2="${cutoffs2} na"
      else
        cutoffs2="${cutoffs2} ${coul_cut}"
      fi
      ;;
    *) echo "Error: ${style} not recognized"
      ;;
  esac
done

read -ra pair_style_array <<< "${pair_style}"
read -ra cutoffs1_array <<< "${cutoffs1}"
read -ra cutoffs2_array <<< "${cutoffs2}"
pair_style_result="${pair_style_array[0]}"

if [[ "${pair_style_result}" == "hybrid/overlay" ]]; then
  for i in $(seq 1 "${#cutoffs1_array[@]}"); do
    if [[ "${cutoffs1_array[$(($i-1))]}" == "na" ]]; then
      pair_style_result="${pair_style_result} ${pair_style_array[$i]}"
    elif [[ "${cutoffs2_array[$(($i-1))]}" == "na" ]]; then
      pair_style_result="${pair_style_result} ${pair_style_array[$i]} ${cutoffs1_array[$(($i-1))]}"
    else
      pair_style_result="${pair_style_result} ${pair_style_array[$i]} ${cutoffs1_array[$(($i-1))]} ${cutoffs2_array[$(($i-1))]}"
    fi
  done
fi

echo "pair_style ${pair_style_result}" >> "${LAMMPS_INPUT_FILE}"
input=$(mktemp 2>/dev/null) || input=/tmp/input$$
trap 'rm -f ${input}' 0 1 2 3 15

dialog --title "Enter the pair_coeff commands (one per line)" --editbox ${input} 10 80 2>>${LAMMPS_INPUT_FILE}

rcut=0
for num in "${cutoffs1_array[@]} ${cutoffs2_array[@]}"; do
  if [[ "${num}" != "na" ]] && [[ "${num}" -gt "${rcut}" ]]; then
    rcut=${num}
  fi
done

if [[ "${rcut}" -eq 0 ]]; then
  rcut=$((4*${a0})) # the 4 is arbitrary
fi

for i in $(seq 1 ${n_types}); do
  for j in $(seq 1 ${n_types}); do
    continue
    echo "pair_write ${i} ${j} 1000 r 0.01 ${rcut} ${i}_${j}_pair_energy.txt \"${i}-${j} Interaction\" ${charges[${i}]} ${charges[$j]}"
  done
done
