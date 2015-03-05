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

cspad_psana_shape = (4, 8, 185, 388)
cspad_geom_shape  = (1480, 1552)

def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='CsPadMaskMaker - mask making, but with a mouse!')
    parser.add_argument('cspad_fnam', type=str, help="filename for the hdf5 cspad image file.")
    parser.add_argument('h5path', type=str, help="hdf5 path for the 2D cspad data.")
    return parser.parse_args()

def unbonded_pixels():
    def ijkl_to_ss_fs(cspad_ijkl):
        """ 
        0: 388        388: 2 * 388  2*388: 3*388  3*388: 4*388
        (0, 0, :, :)  (1, 0, :, :)  (2, 0, :, :)  (3, 0, :, :)
        (0, 1, :, :)  (1, 1, :, :)  (2, 1, :, :)  (3, 1, :, :)
        (0, 2, :, :)  (1, 2, :, :)  (2, 2, :, :)  (3, 2, :, :)
        ...           ...           ...           ...
        (0, 7, :, :)  (1, 7, :, :)  (2, 7, :, :)  (3, 7, :, :)
        """
        if cspad_ijkl.shape != cspad_psana_shape :
            raise ValueError('cspad input is not the required shape:' + str(cspad_psana_shape) )

        cspad_ij = np.zeros(cspad_geom_shape, dtype=cspad_ijkl.dtype)
        for i in range(4):
            cspad_ij[:, i * cspad_psana_shape[3]: (i+1) * cspad_psana_shape[3]] = cspad_ijkl[i].reshape((cspad_psana_shape[1] * cspad_psana_shape[2], cspad_psana_shape[3]))

        return cspad_ij

    mask = np.ones(cspad_psana_shape)

    for q in range(cspad_psana_shape[0]):
        for p in range(cspad_psana_shape[1]):
            for a in range(2):
                for i in range(19):
                    mask[q, p, i * 10, i * 10] = 0
                    mask[q, p, i * 10, i * 10 + cspad_psana_shape[-1]/2] = 0

    mask_slab = ijkl_to_ss_fs(mask)

    import scipy.signal
    mask_pad = scipy.signal.convolve(1 - mask_slab.astype(np.float), np.ones((3, 3), dtype=np.float), mode = 'same') < 1
    return mask_pad

def asic_edges():
    mask_edges = np.ones(cspad_geom_shape)
    mask_edges[:: 185, :] = 0
    mask_edges[:, :: 194] = 0

    mask_edges_pad = scipy.signal.convolve(1 - mask_edges.astype(np.float), np.ones((8, 8), dtype=np.float), mode = 'same') < 1



if __name__ == '__main__':
    args = parse_cmdline_args()

    # load the image
    cspad = h5py.File(args.cspad_fnam, 'r')[args.h5path].value

    mask = np.ones_like(cspad, dtype=np.bool)

    ## Always start by initializing Qt (only once per application)
    app = QtGui.QApplication([])

    ## Define a top-level widget to hold everything
    w = QtGui.QWidget()


    def update_mask_unbonded(state, mask):
        trans = np.fliplr(cspad.T)
        print state
        if state > 0 :
            print 'applying unbonded pixels'
            mask *= unbonded_pixels()
              
            # 
            cspad_max  = cspad.max()
            trans_mask = np.fliplr(mask.T)
            # convert to RGB
            # Set masked pixels to R
            display_data = np.zeros((trans.shape[0], trans.shape[1], 3), dtype = cspad.dtype)
            print display_data, trans.shape, trans_mask.shape
            display_data[:, :, 0] = trans * trans_mask
            display_data[:, :, 1] = trans * trans_mask
            display_data[:, :, 2] = trans + (cspad_max - trans) * ~trans_mask
            
            plot.setImage(display_data)
        else :
            print 'showing the cspad'
            plot.setImage(trans)

    ## Create some widgets to be placed inside
    unbonded_checkbox = QtGui.QCheckBox('unbonded pixels')
    unbonded_checkbox.stateChanged.connect( lambda x: update_mask_unbonded(x, mask) )

    text = QtGui.QLineEdit('enter text')
    listw = QtGui.QListWidget()
    
    plot = pg.ImageView()

    # mouse hover ij value label
    ij_label = QtGui.QLabel()
    disp = 'ss fs {0:5} {1:5}   value {2:6}'.format('-', '-', '-')
    ij_label.setText(disp)
    def mouseMoved(pos):
        img = plot.getImageItem()
        ij = [cspad.shape[0] - 1 - int(img.mapFromScene(pos).y()), int(img.mapFromScene(pos).x())] # ss, fs
        if (0 <= ij[0] < cspad.shape[0]) and (0 <= ij[1] < cspad.shape[1]):
            ij_label.setText('ss fs value: ' + str(ij[0]).rjust(5) + str(ij[1]).rjust(5) + str(cspad[ij[0], ij[1]]).rjust(8) )
    plot.scene.sigMouseMoved.connect(mouseMoved)

    ## Create a grid layout to manage the widgets size and position
    layout = QtGui.QGridLayout()
    w.setLayout(layout)

    ## Add widgets to the layout in their proper positions
    layout.addWidget(ij_label, 0, 0)   # upper-left
    layout.addWidget(unbonded_checkbox, 1, 0)       # text edit goes in middle-left
    layout.addWidget(listw, 2, 0)      # list widget goes in bottom-left
    layout.addWidget(plot, 0, 1, 3, 1) # plot goes on right side, spanning 3 rows
    layout.setColumnStretch(1, 1)
    layout.setColumnMinimumWidth(0, 10)
    
    # display the image
    plot.setImage(np.fliplr(cspad.T))

    ## Display the widget as a new window
    w.show()

    ## Start the Qt event loop
    #app.exec_()
