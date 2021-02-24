#! /usr/bin/env python3

import os, glob, sys #, re
import regex as re
from colorama import Fore, Back, Style
import argparse

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] [options]', description = "Validates simulation output against input and the current directory")
parser.add_argument('infile', help = "The LAMMPS command file (e.g. *.in or in.*)")
parser.add_argument('--logname', default = "log.lammps", help = "The name of the LAMMPS log file. Default: log.lammps")
parser.add_argument('--subfile', nargs = '+', default = "lmp.slurm", help = "The name(s) of the submission script used. Default: lmp.slurm")
args = parser.parse_args()

checkmark = Fore.GREEN + Style.BRIGHT + u'\u2713' + Style.RESET_ALL
fancy_x = Fore.RED + Style.BRIGHT + u'\u2717' + Style.RESET_ALL
fancy_question = Fore.YELLOW + Style.BRIGHT + '?' + Style.RESET_ALL

DIR_REGEX = re.compile(r"(?P<rotation_axis>1[01][01])/(?P<potential>adp|meam|ternary_eam)/((?P<concentration>[0-9]{1,2}(.[0-9]{1})?)at%(?P<solute>[A-Z]([a-z])?)?/)*(?P<misorientation>(?P<angle>[0-9][0-9])degree|sigma7)/T(?P<temperature>[1-9][0-9][05]0)/(?P<radius>large|medium|small|tiny)_r")
DIR_TEMPERATURE_REGEX = re.compile(r"/T(?P<dir_temperature>[1-9][0-9][05]0)")
SUB_SCRIPT_TEMPERATURE_REGEX = re.compile(r"-var T (?P<sub_script_temperature>([1-9][0-9]+|\$\{(?P<temp_var_name>.*?)\}))")
SUB_SCRIPT_CONCENTRATION_FILE_REGEX = re.compile(r"_(?P<sub_script_concentration>[1-9]*)at%(?P<solute>.{1,2}).dat")
SUB_SCRIPT_CONCENTRATION_GENERATE_REGEX = re.compile(r"U (?P<solute_concentration>0.[0-9][0-9]) 1")
SUB_SCRIPT_FILE_REGEX = re.compile(r"-var data_?file (?P<datafile>(LAMMPS_.*\.dat|\$\{.*\}))")
SUB_SCRIPT_LATTICE_PARAM = re.compile(r"\${RANDOM} (?P<lattice_param>3\.[0-9]{2,3})")
DATA_FILE_REGEX = re.compile(r"LAMMPS_(?P<element>([A-Z]([a-z])?([1-9])*)+)_(?P<axis>[0-9][0-9][0-9])_(?P<angle>[1-9][0-9](\.[0-9]+)?)degree_r(?P<radius>[1-9][0-9]*)A_p(?P<potential>([A-Z][a-z]*)+)?(_c(?P<concentration>0\.[0-9]+)(?P<solute>[A-Z][a-z]))?(?P<extra>_[a-zA-Z0-9_.]*)*\.(dat|lmp)")
LOG_FILE_POTENTIAL_REGEX = re.compile(r"pair_style (?P<potential>meam/c|adp|eam/alloy)")

def print_error(var_not_found, string):
    print(f"{Fore.YELLOW}{var_not_found.capitalize()} not found in {string}{Style.RESET_ALL}")

def validate_structure_file(file):
    dirpath = os.path.dirname(os.path.abspath(file))
    dir_radius_dict = {"large":100, "medium":75, "small":50, "tiny":30}
    dir_potential_dict = {"adp":3.52, "ternary_eam":3.542, "meam":3.463}
    file_potential_dict = {"adp":3.52, "eam/alloy":3.542, "meam/c":3.463}
    attribute_dict = {"rotation_axis": dict(), "angle": dict(), "radius": dict(), "potential": dict(), "concentration": dict()}

    # This is something I can try to use, but each case is slightly different, so i
    # will likely need to refactor the code so this can be used.
    # for key in attribute_dict.keys():
    #     if key == "angle":
    #         continue
    #     try:
    #         attribute_dict[key]["dir"] = DIR_REGEX.search(dirpath).group(key)
    #     except AttributeError:
    #         attribute_dict[key]["dir"] = None
    try:
        attribute_dict["rotation_axis"]["dir"] = int(DIR_REGEX.search(dirpath).group('rotation_axis'))
    except AttributeError:
        attribute_dict["rotation_axis"]["dir"] = None
    try:
        attribute_dict["angle"]["dir"] = float(DIR_REGEX.search(dirpath).group('angle'))
    except (AttributeError, TypeError):
        try:
            attribute_dict["angle"]["dir"] = DIR_REGEX.search(dirpath).group('misorientation')
        except AttributeError:
            attribute_dict["angle"]["dir"] = None
    try:
        dir_radius_key = DIR_REGEX.search(dirpath).group('radius')
        attribute_dict["radius"]["dir"] = dir_radius_dict[dir_radius_key]
    except AttributeError:
        attribute_dict["radius"]["dir"] = None
    try:
        dir_potential = DIR_REGEX.search(dirpath).group('potential')
        attribute_dict["potential"]["dir"] = dir_potential_dict[dir_potential]
    except:
        attribute_dict["potential"]["dir"] = None
    try:
        concentrations = DIR_REGEX.search(dirpath).capturesdict()['concentration']
        solutes = DIR_REGEX.search(dirpath).capturesdict(0['solute'])
        if len(concentrations) == 0:
            attribute_dict["concentration"]["dir"] = None
        else:
            attribute_dict["concentration"]["dir"] = {solutes[i]: float(concentrations[i]) for i in range(len(solutes))}
    except:
        attribute_dict["concentration"]["dir"] = None

    try:
        attribute_dict["rotation_axis"]["file"] = int(re.search(r"_(?P<rotation_axis>1[01][01])_", file).group('rotation_axis'))
    except AttributeError:
        attribute_dict["rotation_axis"]["file"] = None
    try:
        attribute_dict["angle"]["file"] = float(re.search(r"_(?P<angle>[0-9][0-9](\.[0-9][0-9]))?degree_", file).group('angle'))
    except AttributeError:
        attribute_dict["angle"]["file"] = None
    try:
        attribute_dict["radius"]["file"] = int(re.search(r"_r([1-9][0-9]*)A", file).group(1))
    except AttributeError:
        attribute_dict["radius"]["file"] = None
    try:
        attribute_dict["potential"]["file"] = float(SUB_SCRIPT_LATTICE_PARAM.search(file).group('lattice_param'))
    except:
        try:
            with open(args.logname) as f:
                for line in f:
                    if len(LOG_FILE_POTENTIAL_REGEX.findall(line)) >= 1:
                        log_file_potential = LOG_FILE_POTENTIAL_REGEX.search(line).group('potential')
                        attribute_dict["potential"]["file"] = file_potential_dict[log_file_potential]
                        break
        except:
            attribute_dict["potential"]["file"] = None
    try:
        attribute_dict["concentration"]["file"] = float(SUB_SCRIPT_CONCENTRATION_GENERATE_REGEX.search(file).group('solute_concentration'))
    except:
        attribute_dict["concentration"]["file"] = None

    # first element for each key is whether or not we will compare the dir and file values
    # second element if they match
    # third is the printed response
    comparison_dict = {"rotation_axis": [True, False, ""], "angle": [True, False, ""], "radius": [True, False, ""],
                        "potential": [True, False, ""], "concentration": [True, False, ""]}
    angle_is_special = False

    # Check for valid path values
    for key,value in attribute_dict.items():

        if key == "angle" and isinstance(value["dir"], str) and value["dir"] == "sigma7":
            if value["file"] == 38.20:
                comparison_dict[key][0] = True
                angle_is_special = True
            else:
                print_error(key, dirpath)
                comparison_dict[key][0] = False
        elif None in value.values():
            comparison_dict[key][0] = False
            if all(i == None for i in value.values()):
                continue
            if value["dir"] is None:
                print_error(key, dirpath)
            elif value["file"] is None:
                print_error(key, file)

    deleted_keys = []
    for key, value in comparison_dict.items():
        cKey = key.capitalize().replace("_", " ") # use the key to create a clean label
        if value[0]: # whether or not to compare the values
            if key == "angle" and angle_is_special:
                value[1] = True
                value[2] = f"\t{checkmark} {cKey}"
            else:
                value[1] = attribute_dict[key]["dir"] == attribute_dict[key]["file"] # the actual comparison
                if value[1]: # if the comparison if true
                    value[2] = f"\t{checkmark} {cKey}"
                else: # if the comparison is false
                    value[2] = f"\t{fancy_x} {cKey} ({attribute_dict[key]['dir']} != {attribute_dict[key]['file']})" if value[0] else f"\t{fancy_question} {cKey}"
        else:
            deleted_keys.append(key)

    for key in deleted_keys:
        del comparison_dict[key]

    return all([i[1] for i in comparison_dict.values()]), '\n'.join([i[2] for i in comparison_dict.values()])

def validate_submission_file(file):
    var_regex = re.compile("^(?P<var_name>[A-Z]+)=(?P<var_value>.+)$")
    var_dict = dict()
    # validate the temperature
    using_variable_T = False
    temp_var_name = None
    file_temperature = ""
    dirpath = os.path.dirname(os.path.abspath(file))
    if not os.path.exists(dirpath + '/' + file):
        file = "../" + file
    try:
        dir_temperature = float(DIR_TEMPERATURE_REGEX.search(dirpath).group('dir_temperature'))
    except AttributeError:
        dir_temperature = None
    with open(file) as f:
        for line in f:
            # check if a variable is defined, store it and its string in a dict
            res = var_regex.search(line)
            if len(var_regex.findall(line)) >= 1:
                var_dict[res.group('var_name')] = res.group('var_value')

            if line.lstrip().startswith("mpirun"):
                try: # get the temperature from '-var T <value>'
                    file_temperature = SUB_SCRIPT_TEMPERATURE_REGEX.search(line).group('sub_script_temperature')
                    if file_temperature.startswith("$"):
                        file_temperature = float(var_dict[SUB_SCRIPT_TEMPERATURE_REGEX.search(line).group('temp_var_name')])
                    else:
                        file_temperature = float(file_temperature)
                except AttributeError:
                    print(f"{Fore.YELLOW}Temperature not found in submission script{Style.RESET_ALL}")

    if using_variable_T:
        file_temperature = float(var_dict[temp_var_name])

    temperature_match = (dir_temperature == file_temperature)
    if temperature_match:
        temperature_response = f"\t{checkmark} Temperature"
    else:
        temperature_response = f"\t{fancy_x} Temperature ({dir_temperature} != {file_temperature})"

    return temperature_match, temperature_response

# check the submission file
#   Are we using the correct structure file? (e.g. the structure file as above, or a modified file with impurities?)
# check the log file
#   Are we using the correct structure file? (read_data <structure_file>)

if __name__ == "__main__":
    dat_files = glob.glob("LAMMPS*.dat")
    if len(dat_files) == 0:
        dat_files = glob.glob("../LAMMPS*.dat")
        if len(dat_files) == 0:
            structure_match, structure_output = False, f"{Fore.YELLOW}No files starting with the format LAMMPS*.dat found{Style.RESET_ALL}"
        else:
            structure_match, structure_output = validate_structure_file(dat_files[0]) # only taking the first one, which may not be what we want. Can we specify the file as an argument?
    else:
        structure_match, structure_output = validate_structure_file(dat_files[0])

    successful_run = False
    for subfile in args.subfile:
        try:
            temperature_match, temperature_output = validate_submission_file(subfile)
            successful_run = True
        except:
            pass

    if not successful_run:
        print(f"{fancy_x} {os.getcwd()}")
        print("No submission scripts found that match " + ", ".join([i for i in args.subfile]))
        sys.exit(0)
    if structure_match and temperature_match:
        print(f"{checkmark} {os.getcwd()}")
    else:
        print(f"{fancy_x} {os.getcwd()}")
        if not structure_match:
            print(structure_output)
        if not temperature_match:
            print(temperature_output)
