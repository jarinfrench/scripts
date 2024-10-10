#! /usr/bin/env python3.6

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import subprocess # for run
from shutil import which # to check if fd or find is available on the system
import os # for getenv
import csv
import itertools # for zip_longest

def shorten_path(path, end_dirs = 4):
    cwd = path.split("/")
    shorten_individual_dirs = False if os.getenv('PROMPT_DIRTRIM') is not None and int(os.getenv('PROMPT_DIRTRIM')) > 0 else True
    user = os.getenv('USER')
    sysname = os.getenv('SYSNAME')
    if os.getenv('PROMPT_DIRTRIM') is not None:
        num_dirs = int(os.getenv('PROMPT_DIRTRIM'))
        if num_dirs <= 3:
            num_dirs = 4
    else:
        num_dirs = 4

    if len(cwd) < 3:
        result = ""
    else:
        if "/".join(cwd[0:3]) == f"/home/{user}":
            result = "~"
            cwd = cwd[3:] # truncate the beginning part of the path
        elif "/".join(cwd[0:3]) == f"/media/{user}":
            result = "m~"
            if len(cwd) > 3:
                dir_num = cwd[3].lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvqxyz./,_!@#$%^&*()-=+:;'\"<>? ") # strips away all letters and characters from the left, leaving only numbers
                if len(cwd[3].split()) > 1: # More than one word in this directory
                    result += cwd[3].split()[0][0:3].upper()
                else: # just one word
                    for i,char in enumerate(cwd[3]):
                        if i >= 3:
                            break
                        else:
                            result += char.upper()
                result += dir_num
                cwd = cwd[4:] # truncate the beginning part of the path
            else:
                cwd = [""]
        elif len(cwd) > 4 and "/".join(cwd[0:4]) == f"/work/{sysname}/{user}":
            result = "w~"
            cwd = cwd[4:]
        else:
            result = ""

    if shorten_individual_dirs:
        if len("/".join([result] + cwd)) < 50:
            return "/".join([result] + cwd)
        for directory in cwd:
            if len(directory) > 8:
                result += "/{}...".format(directory[0:5])
            else:
                result += "/{}".format(directory)
        return result
    else:
        if len(cwd) > num_dirs:
            cwd = [cwd[0]] + ["..."] + cwd[-end_dirs:]
        return "/".join([result] + cwd)

def is_tool(name):
    return which(name) is not None

root = tk.Tk()
root.title("Area Data Comparison")
root.geometry("1500x800")

class AreaDataComparisonApp(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.adc_frame = tk.Frame(master, width = 1200, height = 600)
        self.adc_frame.grid()
        self._plot_styles = ['bo', 'go', 'ro', 'co', 'mo', 'yo', 'ko',
                             'bv', 'gv', 'rv', 'cv', 'mv', 'yv', 'kv',
                             'b^', 'g^', 'r^', 'c^', 'm^', 'y^', 'k^',
                             'b<', 'g<', 'r<', 'c<', 'm<', 'y<', 'k<',
                             'b>', 'g>', 'r>', 'c>', 'm>', 'y>', 'k>',
                             'b1', 'g1', 'r1', 'c1', 'm1', 'y1', 'k1',
                             'b2', 'g2', 'r2', 'c2', 'm2', 'y2', 'k2',
                             'b3', 'g3', 'r3', 'c3', 'm3', 'y3', 'k3',
                              'b4', 'g4', 'r4', 'c4', 'm4', 'y4', 'k4',
                             'b8', 'g8', 'r8', 'c8', 'm8', 'y8', 'k8',
                             'bs', 'gs', 'rs', 'cs', 'ms', 'ys', 'ks',
                             'bp', 'gp', 'rp', 'cp', 'mp', 'yp', 'kp',
                             'bP', 'gP', 'rP', 'cP', 'mP', 'yP', 'kP',
                             'b*', 'g*', 'r*', 'c*', 'm*', 'y*', 'k*',
                             'bh', 'gh', 'rh', 'ch', 'mh', 'yh', 'kh',
                             'bH', 'gH', 'rH', 'cH', 'mH', 'yH', 'kH',
                             'b+', 'g+', 'r+', 'c+', 'm+', 'y+', 'k+',
                             'bx', 'gx', 'rx', 'cx', 'mx', 'yx', 'kx',
                             'bX', 'gX', 'rX', 'cX', 'mX', 'yX', 'kX',
                             'bD', 'gD', 'rD', 'cD', 'mD', 'yD', 'kD',
                             'bd', 'gd', 'rd', 'cd', 'md', 'yd', 'kd',
                             'b|', 'g|', 'r|', 'c|', 'm|', 'y|', 'k|',
                             'b_', 'g_', 'r_', 'c_', 'm_', 'y_', 'k_']

        # variable initialization and defaults
        ## Boolean show_legend keeps track of whether the legend is shown or not
        self.show_legend = tk.IntVar()
        self.show_legend.set(1)
        self.rootdir = tk.StringVar()
        self.rootdir.set("/media/jarinf/Research2/tmp/home/jarinf/projects/uo2/grain_growth/cylindrical")
        self.file_glob = tk.StringVar()
        self.file_glob.set("area_data*")
        self.max_dir_depth = tk.StringVar()
        self.max_dir_depth.set("7") # this is a sensible default based on my current directory structure
        self._last_max_dir_depth = 7
        self.datafiles_frame = tk.Frame(self.adc_frame, width = 350, height = 200)
        self.datafiles_names = []

        ## figure window
        self.fig = plt.figure(figsize=(8,6), dpi=100)

        # All the widgets
        self.show_legend_check_btn = tk.Checkbutton(self.adc_frame, text = "Show legend", variable = self.show_legend, command = self._toggle_legend)
        self.show_legend_check_btn.deselect() # drootdir_labelefault to on

        self.canvas = FigureCanvasTkAgg(self.fig, self.adc_frame)
        self.canvas.draw()

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.adc_frame, pack_toolbar = False)
        self.toolbar.update()

        self.rootdir_update_button = tk.Button(self.adc_frame, text = "Change root directory", command = self._open)
        self.rootdir_label = tk.Label(self.adc_frame, text = "Current directory: " + shorten_path(str(self.rootdir.get())))

        self.file_glob_entry = tk.Entry(self.adc_frame, textvariable = self.file_glob)
        self.file_glob_label = tk.Label(self.adc_frame, text = "File glob: ")
        self.file_glob.trace_add('write', lambda name, index, mode, var = self.file_glob: self._populate(self))

        self.max_dir_depth_entry = tk.Entry(self.adc_frame, textvariable = self.max_dir_depth, validate = 'key', validatecommand = (self.adc_frame.register(self._validate_as_num), '%P'))
        self.max_dir_depth_label = tk.Label(self.adc_frame, text = "Max depth: ")
        self.max_dir_depth.trace_add('write', lambda name, index, mode, var = self.max_dir_depth: self._populate(self))

        self.datafiles_yscrollbar = tk.Scrollbar(self.datafiles_frame, orient = tk.VERTICAL)
        self.datafiles_listbox = tk.Listbox(self.datafiles_frame, width = 50, yscrollcommand = self.datafiles_yscrollbar.set, selectmode = tk.MULTIPLE)
        self.datafiles_listbox.config(exportselection = False)
        self._populate() # initial population of the list
        self.datafiles_listbox.bind('<<ListboxSelect>>', self._plot)
        self.datafiles_yscrollbar.config(command = self.datafiles_listbox.yview)


        # Put everything on the screen
        self.show_legend_check_btn.grid(row = 0, column = 10, sticky = "NW")
        self.canvas.get_tk_widget().grid(row = 0, column = 0, columnspan = 10, rowspan = 10, sticky = "W", pady = 2)
        self.toolbar.grid(row = 10, column = 2, columnspan = 6, sticky = "SEW")
        self.rootdir_update_button.grid(row = 11, column = 0, columnspan = 3, sticky = "W")
        self.rootdir_label.grid(row = 12, column = 0, columnspan = 10, sticky = "W")
        self.file_glob_label.grid(row = 3, column = 11, sticky = "E")
        self.file_glob_entry.grid(row = 3, column = 12, columnspan = 6, sticky = "W")
        self.max_dir_depth_label.grid(row = 3, column = 18, sticky = "E", padx = 5)
        self.max_dir_depth_entry.grid(row = 3, column = 19, columnspan = 3, sticky = "W")
        self.datafiles_yscrollbar.grid(row = 4, column = 23, rowspan = 6, sticky = "ENS")
        self.datafiles_frame.grid(row = 4, column = 11, columnspan = 12, rowspan = 6)
        self.datafiles_listbox.grid(row = 4, column = 11, columnspan = 12, rowspan = 6)

    def _plot(self, event):
        self.fig.clf()
        self.fig.add_subplot(111)
        self.ax = plt.gca()
        data = []
        x_max = -1.0e20
        x_min = 1.0e20
        y_max = -1.0e-20
        if len(self.datafiles_listbox.curselection()) == 0:
            self.ax.clear()
            self.ax.figure.canvas.draw()
            return
        paths = [os.path.splitext(self.datafiles_names[i])[0] for i in self.datafiles_listbox.curselection()]
        labels = [[] for i in range(len(paths))]
        title_vals = []
        for items in itertools.zip_longest(*([i.split('/') for i in paths])):
            if len(set(items)) == 1: # only one unique value
                title_vals.append(items[0])
            else: # multiple values, so include it in the label
                for i in range(len(items)):
                    if items[i] is not None:
                        labels[i].append(items[i])

        title = ' '.join([i for i in paths[0].split('/') if i in title_vals and not i == "area_data"])
        for idx, p in enumerate(paths):
            unique = [i for i in p.split('/') if not i in title_vals]
            labels[idx] = '_'.join(unique)

        if len(labels) == 1:
            labels[0] = "area_data.txt"

        for idx,i in enumerate(self.datafiles_listbox.curselection()):
            x_data = []
            y_data = []
            fin = open(self.rootdir.get() + "/" + self.datafiles_names[i])
            try:
                reader = csv.reader(fin, delimiter = ',')
                for j, line in enumerate(reader):
                    if not (line) or line[0].startswith('#'):
                        continue
                    else:
                        x_data.append(float(line[0])) # this is specific to area_data.txt files
                        y_data.append(float(line[2]))
            except ValueError:
                reader = csv.reader(fin, delimiter = ' ')

                for j,line in enumerate(reader):
                    if not (line) or line[0].startswith("#"):
                        continue
                    else:
                        x_data.append(float(line[0])) # this is specific to area_data.txt files
                        y_data.append(float(line[2]))
            xy1 = sorted(zip(x_data,y_data))
            x_data = [x for x,y in xy1]
            y_data = [y for x,y in xy1]
            x_max = max(max(x_data), x_max)
            x_min = min(min(x_data), x_min)
            y_max = max(max(y_data), y_max)
            plt.plot(x_data, y_data, self._plot_styles[idx], label = labels[idx], markersize = 0.5)

        self.ax.set_xlabel("Time (ps)")
        self.ax.set_ylabel("Area")
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(0, y_max * 1.05)
        self.ax.set_title(title + "\nGrain Area vs Time")
        self.canvas.draw()
        self.ax.legend(loc='best', markerscale = 10)
        self._toggle_legend()

    def _toggle_legend(self):
        h, _ = self.ax.get_legend_handles_labels()
        if len(h) == 0:
            return

        if self.show_legend.get() == 1:
            self.ax.legend(loc='best', markerscale = 10).set_visible(True)
        else:
            self.ax.legend(loc='best', markerscale = 10).set_visible(False)
        self.canvas.draw()

    def _open(self):
        self.rootdir.set(tk.filedialog.askdirectory(initialdir = self.rootdir.get(), title = "Choose a directory", mustexist = True))
        self.rootdir_label['text'] = "Current directory: " + shorten_path(str(self.rootdir.get()))

    # Adapted from https://stackoverflow.com/a/4140988
    def _validate_as_num(self, value):
        if value == '' or value.isdigit():
            if value.isdigit():
                self._last_max_dir_depth = value
            return True
        else:
            self.adc_frame.bell()
            return False

    def _populate(self, _ = "None"):
        self.datafiles_listbox.delete(0,tk.END) # clears the current listbox
        self.datafiles_names = [] # clear the existing file names
        glob = self.file_glob.get()
        if is_tool("fd"):
            if "/" in self.file_glob.get():
                glob = glob.replace("*", "**")
                output = subprocess.run([which("fd"), "-pg", "**/" + glob, "--max-depth", self.max_dir_depth.get(), "-t", "f"], cwd = self.rootdir.get(), stdout = subprocess.PIPE)
            else:
                output = subprocess.run([which("fd"), glob, "--max-depth", self.max_dir_depth.get(), "--glob", "-t", "f"], cwd = self.rootdir.get(), stdout = subprocess.PIPE)
        else:
            if "/" in self.file_glob.get():
                output = subprocess.run([which("find"), ".", "-maxdepth", self.max_dir_depth.get(), "-path", "*" + glob], cwd = self.rootdir.get(), stdout = subprocess.PIPE)
            else:
                output = subprocess.run([which("find"), ".", "-maxdepth", self.max_dir_depth.get(), "-name", glob], cwd = self.rootdir.get(), stdout = subprocess.PIPE)
        files_list = output.stdout.decode("utf-8").split("\n")
        files_list.sort()
        for file in files_list:
            if file:
                self.datafiles_names.append(file)
                self.datafiles_listbox.insert(tk.END, shorten_path(file,4))



def update_plot(event):
    return None

app = AreaDataComparisonApp(root)

root.mainloop()
