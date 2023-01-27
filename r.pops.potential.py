#!/usr/bin/env python3
############################################################################
#
# MODULE:       r.pops.potential
# AUTHOR(S):    Anna Petrasova
# PURPOSE:      Calculating infestation potential
# COPYRIGHT:    (C) 2023 by Anna Petrasova, and the GRASS Development Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
############################################################################


# %module
# % description: Computes infestation potential as a function of susceptible and infected host.
# % keyword: raster
# % keyword: disease
# % keyword: pest
# % keyword: parallel
# %end
# %option G_OPT_R_INPUT
# % key: host
# % label: Input host raster map
# % description: Number of hosts per cell.
# %end
# %option G_OPT_R_INPUT
# % key: infected
# % label: Number of infected hosts
# % description: Number of infected hosts per cell
# %end
# %option G_OPT_R_INPUT
# % key: weather
# % required: no
# % description: Raster of average weather coefficients
# %end
# %option G_OPT_R_OUTPUT
# % key: infestation_potential
# % description: Infestation potential
# %end
# %option G_OPT_R_OUTPUT
# % key: infestation_range
# % required: yes
# % description: Infestation range
# %end
# %option
# % key: reproductive_rate
# % type: double
# % required: yes
# % label: Number of spores or pest units produced by a single host
# %end
# %option
# % key: natural_distance
# % type: double
# % required: yes
# % label: Distance parameter for dispersal
# %end
# %option G_OPT_M_NPROCS
# %end

import atexit
import numpy as np
from math import ceil, sqrt

import grass.script as gs

TMPFILE = None
TMP = []


def cleanup():
    if TMP:
        gs.run_command("g.remove", flags="f", type=["raster"], name=TMP, quiet=True)
    if TMPFILE:
        gs.try_remove(TMPFILE)


def main():
    options, flags = gs.parser()
    natural_distance = float(options["natural_distance"])
    reproductive_rate = float(options["reproductive_rate"])
    infestation_range = options["infestation_range"]
    weather = options["weather"] if options["weather"] else 1

    gs.mapcalc(
        f"{infestation_range} = {natural_distance} * sqrt({options['infected']} * {weather} * {reproductive_rate} * 52 - 1)",
        quiet=True,
    )
    max_distance = gs.raster_info(infestation_range)["max"]
    gs.message(f"Maximum infestation range is {max_distance}")
    matrix = distance_matrix(max_distance, natural_distance)

    path = gs.tempfile()
    global TMPFILE
    TMPFILE = path
    with open(path, "w") as f:
        f.write(write_filter(matrix))

    rmfilter_in = gs.append_node_pid("rmfilter_in")
    TMP.append(rmfilter_in)
    rmfilter_out = gs.append_node_pid("rmfilter_out")
    TMP.append(rmfilter_out)

    gs.mapcalc(
        f"{rmfilter_in} = abs({options['host']} - {options['infected']})", quiet=True
    )
    gs.run_command(
        "r.mfilter",
        input=rmfilter_in,
        output=rmfilter_out,
        filter=path,
        nprocs=int(options["nprocs"]),
        quiet=True,
    )
    gs.mapcalc(
        f"{options['infestation_potential']} = {options['infected']} * {weather} * {rmfilter_out}",
        quiet=True,
    )


def distance_matrix(max_distance, natural_distance):
    region = gs.region()
    size = ceil(max_distance / min(region["ewres"], region["nsres"]))
    resolution = (region["ewres"] + region["nsres"]) / 2
    matrix_size = 2 * size + 1
    matrix = np.zeros((matrix_size, matrix_size))
    center = size
    for i in range(matrix_size):
        for j in range(matrix_size):
            if i != center and j != center:
                dist = (
                    sqrt((i - center) * (i - center) + (j - center) * (j - center))
                    * resolution
                )
                if dist <= max_distance:
                    matrix[i, j] = 1 / (
                        dist * dist + natural_distance * natural_distance
                    )
    return matrix


def write_filter(matrix):
    filter_text = ["TITLE infestation potential"]
    filter_text.append("MATRIX %s" % matrix.shape[0])
    for i in range(matrix.shape[0]):
        line = ""
        for j in range(matrix.shape[0]):
            line += str(matrix[i, j])
            line += " "
        filter_text.append(line)
    filter_text.append("DIVISOR 1")
    filter_text.append("TYPE P")

    return "\n".join(filter_text)


if __name__ == "__main__":
    atexit.register(cleanup)
    main()
