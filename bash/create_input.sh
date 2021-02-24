#! /bin/bash

echo "Specify the type of input file to create."
echo "Note that they are created based on files in the current directory!"
select input_type in find_grains rotate_and_remove calculate_displacement calculate_grain_area generate_impurities grains2slice calculate_plane_density; do
  case ${input_type} in
    find_grains)
      # input data for find_grains:
      # number of files, additional rotation, number of atom types, neighbor cutoff distance, orientation parameter cutoff, lattice parameter
      read -rp "Enter the rotation: " rotation
      read -rp "Enter the number of atom types: " ntypes
      read -rp "Enter the neighbor cutoff distance: " rcut
      read -rp "Enter the orientation parameter cutoff: " fi_cut
      read -rp "Enter the 0K lattice parameter: " a0
      read -rp "Enter the rotation axis: " axis
      read -rp "Enter the crystal structure (bcc|fcc): " crys_struct
      if ! [[ "${crys_struct}" =~ [fb]cc ]]; then
        echo "Crystal structure needs to be either bcc or fcc"
        exit
      fi

      dest_file=$(python3 -c "from myModules import verify_new_file; print(verify_new_file('find_grains_input.txt'))")
      # echo "Running command: 'echo "$(ls -v *0.dump | wc -l) ${rotation} ${ntypes} ${rcut} ${fi_cut} ${a0} ${crys_struct}" > "${dest_file}"'"
      echo "$(find . -maxdepth 1 -name "*0.dump" | wc -l) ${rotation} ${ntypes} ${rcut} ${fi_cut} ${a0} ${crys_struct}" > "${dest_file}"
      case ${axis} in
        100)
          echo "1 0 0" >> "${dest_file}"
          echo "0 1 0" >> "${dest_file}"
          echo "0 0 1" >> "${dest_file}"
          ;;
        110)
          echo "0 0 1" >> "${dest_file}"
          echo "1 -1 0" >> "${dest_file}"
          echo "1 1 0" >> "${dest_file}"
          ;;
        111)
          echo "1 -1 0" >> "${dest_file}"
          echo "1 1 -2" >> "${dest_file}"
          echo "1 1 1" >> "${dest_file}"
          ;;
        112)
          echo "1 1 1" >> "${dest_file}"
          echo "1 -1 0" >> "${dest_file}"
          echo "1 1 -2" >> "${dest_file}"
          ;;
        *)
          read -rp "Please enter the orientation of the x direction (x1 x2 x3): " x1 x2 x3
          read -rp "Please enter the orientation of the y direction (y1 y2 y3): " y1 y2 y3
          read -rp "Please enter the orientation of the z direction (z1 z2 z3): " z1 z2 z3

          echo "${x1} ${x2} ${x3}" >> "${dest_file}"
          echo "${y1} ${y2} ${y3}" >> "${dest_file}"
          echo "${z1} ${z2} ${z3}" >> "${dest_file}"
          ;;
      esac
      ls -v ./*0.dump >> "${dest_file}"
      echo "${dest_file} was created for the find_grains script on $(date)" >> input_file_readme
      break
      ;;
    calculate_displacement)
    # input data for calculate_displacement
      dest_file=$(python3 -c "from myModules import verify_new_file; print(verify_new_file('displacement_input.txt'))")
      ls -v ./*0.dump > "${dest_file}"
      echo "${dest_file} was created for the calculate_displacement script on $(date)" >> input_file_readme
      break
      ;;
    calculate_grain_area)
      # input data for calculate_grain_area
      read -rp "Enter the temperature: " T
      read -rp "Enter the solute concentration in atomic percent: " conc
      read -rp "Enter the 0K lattice parameter: " a0
      read -rp "Enter crystal structure (bcc|fcc): " crys_struct
      echo "Enter the number of the LAMMPS input data file from the list"
      select datafile in *.dat; do break; done
      Lz=$(head -n 7 "${datafile}" | tail -n 1 | awk '{print $2}')
      dest_file=$(python3 -c "from myModules import verify_new_file; print(verify_new_file('grain_area_input.txt'))")
      echo "data.txt ${T} ${conc} ${Lz} ${a0} ${crys_struct}" > "${dest_file}"
      echo "${dest_file} was created for the calculate_grain_area script on $(date)" >> input_file_readme
      break
      ;;
    calculate_plane_density)
      # Input data for calculate_plane_density
      read -rp "Enter the filename containing the list of atoms in the unit cell: " datafile
      if [ ! -f "${datafile}" ]; then
        read -rp "Create the datafile (y|n)? " confirm
        case ${confirm} in
          y|Y)
          echo "Assuming the lattice parameters in each direction are 1.0 Angstroms"
          select struct in fcc bcc sc other; do
            case ${struct} in
              fcc)
              echo "0.0 0.0 0.0" > "${datafile}"
              echo "0.0 0.5 0.5" >> "${datafile}"
              echo "0.5 0.0 0.5" >> "${datafile}"
              echo "0.5 0.5 0.0" >> "${datafile}"
              break
              ;;
              bcc)
              echo "0.0 0.0 0.0" > "${datafile}"
              echo "0.5 0.5 0.5" >> "${datafile}"
              break
              ;;
              sc)
              echo "0.0 0.0 0.0" > "${datafile}"
              break
              ;;
              *)
              read -rp "Enter the number of atoms in the unit cell" num_atoms
              for i in $(seq 1 "${num_atoms}"); do
                read -rp "Enter the coordinates of atom ${i}: " x y z
                echo "$x $y $z" >> "${datafile}"
              done
              break
              ;;
            esac
            echo "${datafile} was created for the calculate_planar_density script on $(date)" >> input_file_readme
          done
          ;;
          n|N)
          echo "Did not create ${datafile}"
          ;;
          *)
          echo "Unrecognized option. Please enter y|n next time."
          exit
          ;;
        esac
      fi
      read -rp "Enter the a lattice constant: " a
      read -rp "Enter the b lattice constant: " b
      read -rp "Enter the c lattice constant: " c
      read -rp "Enter the alpha angle: " alpha
      read -rp "Enter the beta angle: " beta
      read -rp "Enter the gamma angle: " gamma

      dest_file=$(python3 -c "from myModules import verify_new_file; print(verify_new_file('calculate_planar_density_input.txt'))")
      echo "${datafile} $a $b $c ${alpha} ${beta} ${gamma}" > "${dest_file}"
      echo "${dest_file} was created for the calculate_planar_density script on $(date)" >> input_file_readme
      ;;
    rotate_and_remove|generate_impurities|grains2slice)
      echo "Not implemented yet"
      exit
      ;;
    *)
      echo "No input file created: unrecognized value"
      exit
      ;;
  esac
done
