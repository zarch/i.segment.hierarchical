#!/usr/bin/env python
# -- coding: utf-8 --
#
############################################################################
#
# MODULE:	    i.segment.hierarchical
#
# AUTHOR(S):   Pietro Zambelli (University of Trento)
#
# COPYRIGHT:	(C) 2013 by the GRASS Development Team
#
#		This program is free software under the GNU General Public
#		License (>=v2). Read the file COPYING that comes with GRASS
#		for details.
#
#############################################################################

#%Module
#%  description: Hierarchical segmentation
#%  keywords: imagery
#%  keywords: segment
#%  overwrite: yes
#%End
#%option G_OPT_I_GROUP
#%  key: group
#%  description: Name of input imagery group
#%  required: yes
#%end
#%option
#%  key: thresholds
#%  type: double
#%  multiple: yes
#%  description: Segment thresholds
#%  required: yes
#%  answer: 0.02,0.05
#%end
#%option G_OPT_R_OUTPUT
#%  key: output
#%  description: Name of output sement raster map
#%  required: yes
#%end
#%option
#%  key: outputs_prefix
#%  type: string
#%  description: Name for output raster maps from segment
#%  required: no
#%  answer: seg__%.2f
#%end
#%option
#%  key: method
#%  type: string
#%  description: Segmentation method
#%  required: no
#%  answer: region_growing
#%  guisection: Segment
#%end
#%option
#%  key: similarity
#%  type: string
#%  description: Similarity calculation method
#%  required: no
#%  answer: euclidean
#%  guisection: Segment
#%end
#%option
#%  key: minsizes
#%  type: integer
#%  description: Minimum number of cells in a segment
#%  multiple: yes
#%  required: no
#%  guisection: Segment
#%end
#%option
#%  key: memory
#%  type: integer
#%  description: Memory in MB
#%  required: no
#%  answer: 300
#%  guisection: Segment
#%end
#%option
#%  key: iterations
#%  type: integer
#%  description: Maximum number of iterations
#%  required: no
#%  answer: 20
#%  guisection: Segment
#%end
#%option G_OPT_R_INPUT
#%  key: seeds
#%  description: Name for input raster map with starting seeds
#%  required: no
#%  guisection: Segment
#%end
#%option G_OPT_R_INPUT
#%  key: bounds
#%  description: Name of input bounding/constraining raster map
#%  required: no
#%  guisection: Segment
#%end
#%option
#%  key: width
#%  description: Tile width in pixels
#%  type: integer
#%  required: no
#%  guisection: Grid
#%end
#%option
#%  key: height
#%  description: Tile height in pixels
#%  type: integer
#%  required: no
#%  guisection: Grid
#%end
#%option
#%  key: overlap
#%  description: Tile overlap in pixels
#%  type: integer
#%  required: no
#%  answer: 0
#%  guisection: Grid
#%end
#%option
#%  key: processes
#%  description: Number of concurrent processes
#%  type: integer
#%  required: no
#%  guisection: Grid
#%end

# move=None, log=False, start_row=0, start_col=0, out_prefix=''

from __future__ import print_function
import multiprocessing as mltp
import time
from grass.script.core import parser
from grass.pygrass.modules import Module
from grass.pygrass.modules.grid import GridModule
from grass.pygrass.modules.grid.split import split_region_tiles
from grass.pygrass.modules.grid.patch import get_start_end_index
from grass.pygrass.raster import RasterRow

DEBUG = False

def rpatch_row(rast, rasts, bboxes, max_rasts):
    """Patch a row of bound boxes."""
    sei = get_start_end_index(bboxes)
    # instantiate two buffer
    buff = rasts[0][0]
    rbuff = rasts[0][0]
    r_start, r_end, c_start, c_end = sei[0]
    for row in xrange(r_start, r_end):
        for col, ras in enumerate(rasts):
            r_start, r_end, c_start, c_end = sei[col]
            buff = ras.get_row(row, buff)
            rbuff[c_start:c_end] = buff[c_start:c_end] + max_rasts[col]
        rast.put_row(rbuff)


def rpatch_map(raster, mapset, mset_str, bbox_list, overwrite=False,
               start_row=0, start_col=0, prefix=''):
    """Patch raster using a bounding box list to trim the raster."""
    #import ipdb; ipdb.set_trace()
    # Instantiate the RasterRow input objects
    rast = RasterRow(prefix + raster, mapset)
    with RasterRow(name=raster, mapset=mset_str % (0, 0), mode='r') as rtype:
        rast.open('w', mtype=rtype.mtype, overwrite=overwrite)
    rasts = []
    mrast = 0
    for row, rbbox in enumerate(bbox_list):
        rrasts = []
        max_rasts = []
        for col in range(len(rbbox)):
            rrasts.append(RasterRow(name=raster,
                                    mapset=mset_str % (start_row + row,
                                                       start_col + col)))
            rrasts[-1].open('r')
            mrast += rrasts[-1].info.max + 1
            max_rasts.append(mrast)
        rasts.append(rrasts)
        rpatch_row(rast, rrasts, rbbox, max_rasts)

    for rrast in rasts:
        for rast_ in rrast:
            rast_.close()
    rast.close()


class SegModule(GridModule):
    def patch(self):
        """Patch the final results."""
        # patch all the outputs
        bboxes = split_region_tiles(width=self.width, height=self.height)
        inputs = self.module.inputs
        rpatch_map(inputs.outputs_prefix % inputs.thresholds[-1],
                   self.mset.name, self.msetstr, bboxes,
                   self.module.flags.overwrite,
                   self.start_row, self.start_col, self.out_prefix)


def segment(thresholds, minsizes, output='seg__%.2f', **opts):
    iseg = Module('i.segment')
    seeds = None
    for thr, msize in zip(thresholds, minsizes):
        opts['threshold'] = thr
        opts['minsize'] = msize
        opts['seeds'] = seeds
        opts['flags'] = flags
        opts['output'] = output % thr
        st = time.time()
        iseg(**opts)
        print("%s, required: %.2fs" % (opts['output'], time.time() - st))
        seeds = opts['output']  # update seeds


if __name__ == "__main__":
    import os
    GISRC = os.getenv('GISRC')
    opts, flags = parser()
    width = opts.pop('width')
    height = opts.pop('height')
    overlap = opts.pop('overlap')
    processes = opts.pop('processes')
    processes = int(processes) if processes else mltp.cpu_count()
    memory = int(opts['memory'])
    thrs = [float(thr) for thr in opts['thresholds'].split(',')]
    print(repr(opts['minsizes']))
    if opts['minsizes']:
        minsizes = [int(m) for m in opts['minsizes'].split(',')]
        if len(minsizes) != len(thrs):
            minsizes = [int(minsizes[0]), ] * len(thrs)
    else:
        minsizes = [1, ] * len(thrs)

    # define new cleaned parameters
    opts['thresholds'] = thrs
    opts['minsizes'] = minsizes
    opts['iterations'] = int(opts['iterations'])
    opts['memory'] = memory / processes
    if width and height:
        seg = SegModule('i.segment.hierarchical',
                        width=int(width), height=int(height),
                        overlap=int(overlap),
                        processes=processes,
                        debug=DEBUG, **opts)
        #import ipdb; ipdb.set_trace()
        seg.run()
        os.environ['GISRC'] = GISRC
        print("Start running segment for the last time in the whole region")
        st = time.time()
        #import ipdb; ipdb.set_trace()
        iseg = Module('i.segment')
        iseg(group=opts['group'],
             output=opts['output'],
             threshold=opts['thresholds'][-1],
             method=opts['method'],
             similarity=opts['similarity'],
             minsize=opts['minsizes'][-1],
             memory=memory,
             iterations=3,
             seeds=opts['outputs_prefix'] % opts['thresholds'][-1])
        print("%s, required: %.2fs" % (opts['output'], time.time() - st))
    else:
        opts.pop('output')
        segment(opts.pop('thresholds'), opts.pop('minsizes'),
                output=opts.pop('outputs_prefix'), **opts)
