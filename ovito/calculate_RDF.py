#! /usr/bin/env python3

from ovito.io import import_file, export_file
import ovito.modifiers as om
import numpy as np

a0 = 3.542
nn17 = 2 * np.sqrt(3) # 17th nearest neighbor distance in BCC. Arbitrary, but gets us to ~12 angstroms
pipeline = import_file("*.dump")
pipeline.modifiers.append(om.CoordinationAnalysisModifier(cutoff = a0 * nn17, number_of_bins = 200, partial = True))
#pipeline.modifiers.append(om.CoordinationAnalysisModifier(cutoff = a0 * nn17, number_of_bins = 200))

data = pipeline.compute()
nframes = pipeline.source.num_frames
rdf_table = data.tables['coordination-rdf']

export_file(rdf_table, 'rdf_initial.txt', 'txt/table', frame = 0)
export_file(rdf_table, 'rdf_final.txt', 'txt/table', frame = nframes)

#rdf_names = rdf_table.y.component_names
#overall_rdf = data.tables['coordination-rdf.2']
