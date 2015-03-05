#!/usr/bin/env python
"""
 Load a cspad image from the commandline arguement
 Show pixel location from the mouse position
 Show pixel value from the mouse position
 Clicking should add to the mask
 Some sort of ROI?
"""

import argparse
import h5py
from PyQt4 import QtGui
import pyqtgraph as pg
import numpy as np

def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='CsPadMaskMaker - mask making, but with a mouse!')
    parser.add_argument('cspad_fnam', type=str, help="filename for the hdf5 cspad image file.")
    parser.add_argument('h5path', type=str, help="hdf5 path for the 2D cspad data.")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_cmdline_args()

    # load the image
    cspad = h5py.File(args.cspad_fnam, 'r')[args.h5path].value

    ## Always start by initializing Qt (only once per application)
    app = QtGui.QApplication([])

    ## Define a top-level widget to hold everything
    w = QtGui.QWidget()

    ## Create some widgets to be placed inside
    btn = QtGui.QPushButton('press me')
    text = QtGui.QLineEdit('enter text')
    listw = QtGui.QListWidget()

    plot = pg.ImageView()

    ## Create a grid layout to manage the widgets size and position
    layout = QtGui.QGridLayout()
    w.setLayout(layout)

    ## Add widgets to the layout in their proper positions
    layout.addWidget(btn, 0, 0)   # button goes in upper-left
    layout.addWidget(text, 1, 0)   # text edit goes in middle-left
    layout.addWidget(listw, 2, 0)  # list widget goes in bottom-left
    layout.addWidget(plot, 0, 1, 3, 1)  # plot goes on right side, spanning 3 rows
    layout.setColumnStretch(1, 1)
    layout.setColumnMinimumWidth(0, 10)
    
    # display the image
    plot.setImage(np.fliplr(cspad.T))

    ## Display the widget as a new window
    w.show()

    ## Start the Qt event loop
    app.exec_()
