#!/usr/bin/python3

"""Parser for CCP4 files."""

import mrcfile

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import percentileofscore
from numpy import float32

import time
import pandas as pd
import matplotlib.pyplot as plt


class ElectronDensity(object):
    """Parse a CCP4 file and return an ElectronDensityData object."""
    def __init__(self):
        """Create an ElectronDensityData object.
        An ElectronDensityData object is constructed by several parameters that denotes key information usually contained in the header of an CCP4 file, as well as the main data consists of electron density values stored numpy array."""

        # alpha, beta, gamma angles of the cell axes.
        self.cellb = None        

        # number of columns, rows and sections in 3D data array.
        self.column_number = 0
        self.row_number = 0
        self.section_number = 0

        # location of first column, row and section in unit cell.
        self.start_column_number = 0
        self.start_row_number = 0
        self.start_section_number = 0

        # axis corresp to column, row and section (1,2,3 for X,Y,Z)
        self.axis_to_column = 0
        self.axis_to_row = 0
        self.axis_to_section = 0

        # voxel size of unit cell in column, row and section dimension
        self.column_voxel_size = 0
        self.row_voxel_size = 0
        self.section_voxel_size = 0

        # angstrom coordinates for the first column, row and section in unit cell.
        self.start_column_coord = 0
        self.start_row_coord = 0
        self.start_section_coord = 0

        # electron density data in numpy array format, usually appeared in the sequence of section-row-column.
        self.electron_density_data = None

    def fetch_data(self, pdb_id, filepath=None):
        """Take the pdb_id or CCP4 filepath, reads a CCP4 file and return an ElectronDensityData object storing the data in CCP4 file."""
        
        if filepath==None:
            filepath = '../CCP4/{}_2fofc.ccp4'.format(pdb_id)

        try:
            with mrcfile.open(filepath) as mrc_f:
                if mrc_f.data.shape == (mrc_f.header.nz.item(), \
                                        mrc_f.header.ny.item(), \
                                        mrc_f.header.nx.item()):

                    # make sure the data array sequence is section-row-column

                    self.cellb = mrc_f.header.cellb

                    self.column_number = mrc_f.header.nx
                    self.row_number = mrc_f.header.ny
                    self.section_number = mrc_f.header.nz

                    self.start_column_number = mrc_f.header.nxstart
                    self.start_row_number = mrc_f.header.nystart
                    self.start_section_number = mrc_f.header.nzstart

                    self.axis_to_column = mrc_f.header.mapc # 1-x; 2-y; 3-z axis.
                    self.axis_to_row = mrc_f.header.mapr # 1-x; 2-y; 3-z axis.
                    self.axis_to_section = mrc_f.header.maps # 1-x; 2-y; 3-z axis.

                    x_voxel_size = mrc_f.voxel_size.x # voxel size in x axis
                    y_voxel_size = mrc_f.voxel_size.y # voxel size in y axis
                    z_voxel_size = mrc_f.voxel_size.z # voxel size in z axis

                    self.column_voxel_size = [x_voxel_size, y_voxel_size, z_voxel_size]\
                                             [self.axis_to_column-1]
                    self.row_voxel_size = [x_voxel_size, y_voxel_size, z_voxel_size]\
                                          [self.axis_to_row-1]
                    self.section_voxel_size = [x_voxel_size, y_voxel_size, z_voxel_size]\
                                              [self.axis_to_section-1]

                    self.start_column_coord = self.start_column_number * self.column_voxel_size
                    self.start_row_coord = self.start_row_number * self.row_voxel_size
                    self.start_section_coord = self.start_section_number * self.section_voxel_size

                    # arr_sorted = sorted(mrc_f.data.flatten())
                    # data_s = pd.Series(mrc_f.data.flatten())
                    # self.electron_density_data = data_s.apply(lambda x: percentileofscore(arr_sorted, x)).values.reshape(data.shape) # re-normalized standard deviation
                    self.electron_density_data = mrc_f.data # numpy array

                else:
                    print("ERROR! INCONSISTENT DATA SHAPE IN HEADER AND BODY.")
        
        except:
            print("ERROR! NO CCP4 FILE FOUND.")

        return self

    def build_interpolator(self):
        """Interpolation of eletron_density_data onto 3D mesh grid, in order to map actual column, row, and section axis to the absolute electron density value. """
        column_start = self.start_column_coord
        row_start = self.start_row_coord
        section_start = self.start_section_coord

        column_end = column_start + self.column_voxel_size * (self.column_number-1)
        row_end = row_start + self.row_voxel_size * (self.row_number-1)
        section_end = section_start + self.section_voxel_size * (self.section_number-1)

        column_coords = np.linspace(start=column_start, stop=column_end, \
                                    num=self.column_number, endpoint=True)
        row_coords = np.linspace(start=row_start, stop=row_end, \
                                 num=self.row_number, endpoint=True)
        section_coords = np.linspace(start=section_start, stop=section_end, \
                                     num=self.section_number, endpoint=True)

        interpolator = RegularGridInterpolator((section_coords, row_coords, column_coords), \
                                               self.electron_density_data,
                                               bounds_error=False, fill_value=None)
        
        return interpolator

    def axes_conversion(self, xyz_coords):
        """if the cellb of ccp4 file is not (90, 90, 90), we need to convert x0, y0, z0 in orthogonal coordinate system (actual x, y, z in PDB file) to new x1, y1, z1 in non-orghogonal coordinate system."""
        
        alpha, beta, gamma = np.deg2rad(self.cellb.alpha), np.deg2rad(self.cellb.beta), np.deg2rad(self.cellb.gamma)

        old_xyz_coords = np.array(xyz_coords)

        nonortho2ortho_martrix = np.array([[1, 0, 0], [np.cos(gamma), np.sin(gamma), 0], [np.cos(beta), ((np.cos(alpha) - np.cos(beta) * np.cos(gamma))/np.sin(gamma)), np.sqrt((1 - (np.cos(beta))**2 - ((np.cos(alpha) - np.cos(beta) * np.cos(gamma))/np.sin(gamma))**2))]])

        ortho2nonortho = np.linalg.inv(nonortho2ortho_martrix)

        new_xyz_coords = np.dot(old_xyz_coords, ortho2nonortho)
        return new_xyz_coords


    def coord_xyz_to_src(self, xyz_coords):
        """converts x,y,z coordinates (list of sets/tuples) to section,row,column coordinates (list of sets/tuples)."""
        xyz_to_src_map = [self.axis_to_section-1, self.axis_to_row-1, self.axis_to_column-1]

        xyz_coords_array = np.array(xyz_coords)

        new_xyz_coords_array = self.axes_conversion(xyz_coords_array)

        src_coords_array = new_xyz_coords_array[:, xyz_to_src_map]

        return src_coords_array

    def get_density(self, xyz_coords):
        """takes the x,y,z coordinates, return their absolute electron density values as numpy array."""
        interpolator = self.build_interpolator()
        src_coords = self.coord_xyz_to_src(xyz_coords)

        densities = interpolator(src_coords)

        return densities

    def get_maxmin_density(self, atom1_xyz_coord, atom2_xyz_coord):
        
        xyz_coord1 = np.array(atom1_xyz_coord)
        xyz_coord2 = np.array(atom2_xyz_coord)
        line_vector = (xyz_coord2-xyz_coord1)/np.linalg.norm(xyz_coord2-xyz_coord1)
        line_distance = np.linspace(start=0, stop=1, num=100+1, endpoint=True)*np.linalg.norm(xyz_coord2-xyz_coord1)
        line_xyz_coords = (xyz_coord1 + line_vector*line_distance[:, None])

        line_density = self.get_density(line_xyz_coords)

        return line_density.max(), line_density.min()

    def get_line_density(self, atom1_xyz_coord, atom2_xyz_coord, line_start=-2.0, line_end=6.0, relative_ratio=False):
        """get the linear density between any two atoms' xyz coordinates. By default extend the plot by 50% at the start and end of the line connecting the atoms for better visualization"""
        xyz_coord1 = np.array(atom1_xyz_coord)
        xyz_coord2 = np.array(atom2_xyz_coord)

        line_vector = (xyz_coord2-xyz_coord1)/np.linalg.norm(xyz_coord2-xyz_coord1)

        if relative_ratio==False:
            line_distance = np.linspace(start=line_start, stop=line_end, num=100+1, endpoint=True)
        elif relative_ratio==True:
            line_distance = np.linspace(start=line_start, stop=line_end, num=100+1, endpoint=True)*np.linalg.norm(xyz_coord2-xyz_coord1)

        line_xyz_coords = (xyz_coord1 + line_vector*line_distance[:, None])

        line_density = self.get_density(line_xyz_coords)

        return line_distance, line_density

    def get_box_density(self, atom_xyz_coord, boxsize=2, num_interval=10):
        """take atom's coordinates, returns the electron density values of the 2x2x2 Ang^3 box it is in."""
        xyz_coord = np.array(atom_xyz_coord)
        x_coord, y_coord, z_coord = xyz_coord[0], xyz_coord[1], xyz_coord[2]

        boxgrid = np.ones([num_interval, num_interval, num_interval,3])

        linx = np.linspace(start=x_coord-boxsize, stop=x_coord+boxsize, num=num_interval, endpoint=False)
        liny = np.linspace(start=y_coord-boxsize, stop=y_coord+boxsize, num=num_interval, endpoint=False)
        linz = np.linspace(start=z_coord-boxsize, stop=z_coord+boxsize, num=num_interval, endpoint=False)

        xs, ys, zs = np.meshgrid(linx, liny, linz, indexing='ij')

        boxgrid[:,:,:,0] = boxgrid[:,:,:,0]*xs
        boxgrid[:,:,:,1] = boxgrid[:,:,:,1]*ys
        boxgrid[:,:,:,2] = boxgrid[:,:,:,2]*zs

        boxcoords = boxgrid.reshape(1000, 3)
        boxdensities = self.get_density(boxcoords)

        # normalize density data.
        max_den = boxdensities.max()
        min_den = boxdensities.min()
        boxdensities_normed = (boxdensities-min_den) / (max_den-min_den)

        return boxdensities_normed

    def get_mid_density(self, atom1_xyz_coord, atom2_xyz_coord):
        """get the linear density between any two atoms' xyz coordinates. By default extend the plot by 50% at the start and end of the line connecting the atoms for better visualization"""
        xyz_coords_mid = (np.array(atom1_xyz_coord)+np.array(atom2_xyz_coord))/2
        mid_density = self.get_density([xyz_coords_mid, ])[0]

        return mid_density


    def plot_density(self, atom1_xyz_coord, atom2_xyz_coord):
        """plot the linear density between any two atoms' xyz coordinates. By default extend the plot by plot (-2.0 to 6.0 Angstron) at the start and end of the line connecting the atoms for better visualization"""

        line_distance, line_density = self.get_line_density(atom1_xyz_coord, atom2_xyz_coord, line_start=-2.0, line_end=6.0, 
            relative_ratio=False)

        plt.plot(line_distance, line_density)
        plt.xlabel('distance from atom1 to atom 2 ')
        plt.ylabel('electron density from atom1 to atom 2')
        # plt.show()
    # def plot_regularized_density(self, atom1_xyz_coord, atom2_xyz_coord)