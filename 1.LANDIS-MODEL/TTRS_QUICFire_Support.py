'''
Tall Timbers Research Station

'''

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from jinja2 import Template, TemplateError
from scipy import interpolate
from scipy.io import FortranFile

def import_fortran_dat_file(filename,cell_nums):
    (nx,ny,nz) = cell_nums
    arr =  FortranFile(filename,'r','uint32').read_ints('float32').T.reshape(nz,ny,nx)
    return arr

def import_topo_dat_file(filename,cell_nums):
    (nx,ny,nz) = cell_nums
    return FortranFile(filename,'r','uint32').read_ints('float32').T.reshape(ny,nx)
    
def export_fortran_dat_file(arr,filename):
    arr = arr.astype('float32')
    arrfile = FortranFile(filename,'w','uint32')
    arrfile.write_record(arr)

def calculate_QF_vertical_wind_grid(max_topo): 

    c = 1.41/2.0                # 1/sqrt(2)
    e = 2.71828/2.0             # Eulers Number / 2
    zo = e**c                   # Initial height for log profile      
    dist = zo                   # Keeps track of the max height
    zt = zo                     # Keeps track of previous cell's height
    zs = []                     # List storing the cell heights
    zs.append(zo)
    max_topo_adjusted = np.max([max_topo*3.,400.]) # Covering edge scenario for flat topo
    while(dist < max_topo_adjusted):
        zo = 1.125*zt*e**c
        dist+=zo
        zt = zo
        zs.append(zo)
        
    print (zs,len(zs))
    
    return zs
        
def calculate_QF_vertical_fuel_grid(max_fuel_height):  
    # Deprecated, need to move towards using uniform grid
    heights = [1,1,1,1,1.25,1.5,1.75,2.,2.5,3.,3.5,4.,4.5,5.,6.,7.5,9.,12.,15.,15.,15.,15.]
    tot_height = 0
    cnt = 0
    while (np.floor(tot_height) < max_fuel_height):
        tot_height += heights[cnt]
        cnt += 1
    zs = heights[:cnt]
    return zs

def raster_import(filepath):  # NEEDS CLEANING UP
    import gdal
    #Convert raster to np_array
    #https://automating-gis-processes.github.io/2016/Lesson7-read-raster-array.html
    # Open the file:
    raster = gdal.Open(filepath)
    band = raster.GetRasterBand(1)
    
    #Remove null values
    rasterArray = raster.ReadAsArray()
    nodata = band.GetNoDataValue()
    #Remove rows and columns that contain no data
    rasterArray = rasterArray[:,~np.all(rasterArray==nodata, axis=0)] 
    rasterArray = rasterArray[~np.all(rasterArray==nodata, axis=1),:]
    
    ##Convert null to 0
    #rasterArray[np.where(rasterArray==nodata)]=0
    rasterArray = np.flip(rasterArray, axis = 0)
    
    return rasterArray

def interpolate_array(rasterArray, nodata=None,method='nearest'):  # NEEDS CLEANING UP
    # Interpolates from a lower resolution to higher resolution
    x = np.arange(0, rasterArray.shape[1])
    y = np.arange(0, rasterArray.shape[0])
    #mask invalid values
    rasterArray = np.ma.masked_equal(rasterArray, nodata)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~rasterArray.mask]
    y1 = yy[~rasterArray.mask]
    newarr = rasterArray[~rasterArray.mask]
    
    rasterArray = interpolate.griddata((x1, y1), newarr.ravel(),
                              (xx, yy),
                             method=method)
    return rasterArray

def aggregate_array_to_lower_resolution(arr,res_x=2,res_y=2): # NEEDS CLEANING UP
    # Assumes a 'standard' QF dat array (nz,ny,nx)    
    ny,nx  = arr.shape
    arr = arr.reshape((ny//res_y,res_y, nx//res_x,res_x))
    arr = arr.mean(3).mean(1)
    return arr

def remove_shapefile_from_bbox(mask_shp,dom_shp,res_xy=[2,2]):
    '''
    Uses mask shapefile to create a mask in the size of the domain shapefile.
    If you don't have a domain shapefile, feed in the mask_shp.buffer(50) as a good default
   
    The trick for a different routine would be to get the total bounds and build 
    an array of the right 'size' and then feed it in to the larger array. Getting 
    the corner indices will be a little weird.

    '''
    import geopandas as gpd
    from rasterio.features import geometry_mask
    
    mask_shp = gpd.clip(mask_shp, dom_shp)
    x1_min,y1_min,x1_max,y1_max = dom_shp.total_bounds
    x2_min,y2_min,x2_max,y2_max  = mask_shp.total_bounds
    # Defining a bunch of parameters that will be useful later on for cleanness
    msk_xext,msk_yext,dom_xext,dom_yext,x_corner,y_corner = [
        int(x2_max-x2_min),int(y2_max-y2_min),int(x1_max-x1_min),int(y1_max-y1_min),
        int(x2_min-x1_min),int(y2_min-y1_min)]
       
    geoT = (-1.0,0.0,x2_max,
           0.0,-1.0,y2_max)    
    msk = np.fliplr(np.flipud(geometry_mask(mask_shp.geometry, [msk_yext,msk_xext], geoT)))
    arr = np.ones([dom_yext,dom_xext],dtype=int)
    arr[y_corner:y_corner+msk_yext,x_corner:x_corner+msk_xext] *= msk
    arr = arr.reshape(dom_yext//res_xy[1],res_xy[1], dom_xext//res_xy[0],res_xy[0])
    arr = arr.mean(3).mean(1)
    return arr;

class AlbersEqualAreaConic: # Used exclusively to convert Albers  to Lat/Long (Forward) and the Inverse
    

    """
    Implements forward and inverse projection on Albers Equal Area Conic
    Attributes:
        lambda_0 (float): origin of longitude
        a (float): semi-major axis of earth (GRS80)
        e2 (float): eccentricity squared
        e (float): eccentricity
        n (float): stored map constant
        C (float): stored map constant
        rho_0 (float): stored map constant
    """

    def __init__(self, phi_1=29.5, phi_2=45.5, phi_0=23.0, lambda_0=-96.0):
        """
        Creates an instance of AlbersEqualAreaConic, initializes attributes and
        computes map constants
        Note:
            Geographic constants are based on the GRS 1980 ellipsoid
        Args:
            phi_1 (float): first standard parallel
            phi_2 (float): second standard parallel
            phi_0 (float): origin of latitude
            lambda_0 (float): origin of longitude
        """

        # convert map params to radians
        phi_0 = np.radians(phi_0)
        phi_1 = np.radians(phi_1)
        phi_2 = np.radians(phi_2)
        self.lambda_0 = np.radians(lambda_0)

        # GRS 1980 REFERENCE ELLIPSIOD CONSTANTS
        # geographic constants
        self.a = 6378137.0
        # derived geometrical constants
        f = 1.0/298.2572221010042 # flattening
        self.e2 = 2*f - f**2 # eccentricity squared
        self.e = np.sqrt(self.e2) # eccentricity

        # preliminaries
        m_1 = self._m(phi_1)
        m_2 = self._m(phi_2)
        q_0 = self._q(phi_0)
        q_1 = self._q(phi_1)
        q_2 = self._q(phi_2)

        # map constants
        self.n = (m_1**2 - m_2**2)/(q_2 - q_1)
        self.C = m_1**2 + self.n*q_1
        self.rho_0 = self.a*(self.C - self.n*q_0)**0.5/self.n

    def _m(self, phi):
        """Private member
        Convenience method for computing map constants
        """

        return np.cos(phi)/np.sqrt(1 - self.e2*(np.sin(phi))**2)

    def _q(self, phi):
        """Private member
        Another convenience method for computing map constants
        """

        return (1 - self.e2)*(np.sin(phi)/(1 - self.e2*(
            np.sin(phi))**2) - (1.0/(2*self.e))*np.log((1-self.e*np.sin(
            phi))/(1 + self.e*np.sin(phi))))

    def forward(self, lat, lon):
        """
        Performs forward projection from geodetic coordinates to projected
        coordinates
        Args:
            lat (float): latitude
            lon (float): longitude
        Returns:
            (x,y) coordinate projected in Albers Equal Area Conic
        """

        # convert to radians for numpy trig functions
        lat = np.radians(lat)
        lon = np.radians(lon)

        # preliminaries
        q = self._q(lat)
        rho = self.a*(self.C - self.n*q)**0.5/self.n
        theta = self.n*(lon - self.lambda_0)

        # retrieve the projected coordinates
        x = rho*np.sin(theta)
        y = self.rho_0 - rho*np.cos(theta)

        return x,y

    def inverse(self, x, y):
        """
        Performs inverse projection from Albers to geodetic coordinates
        Args:
            x (float): x projected in Albers
            y (float): y projected in Albers
        Returns:
            lat and lon in geodetic coordinates
        """

        # preliminaries
        p = np.sqrt(x*x + (self.rho_0 - y)**2)
        theta = np.arctan2(x, self.rho_0 - y)
        q = (self.C - ((p*p)*(self.n**2))/(self.a**2))/self.n

        # convergence criteria
        epsilon = 1e-6

        # iterate latitude calculation until convergence
        phi = np.sin(q/2)
        next_phi = self._inverse_iteration(phi, q)
        while (np.abs(np.degrees(phi) - np.degrees(next_phi)) > epsilon):
            phi = next_phi
            next_phi = self._inverse_iteration(phi, q)

        return np.degrees(phi), np.degrees(self.lambda_0 + theta/self.n)

    def _inverse_iteration(self, phi, q):
        """Private member
        Formula to iterator until convergence for inverse projection of latitude
        """

        return np.radians(np.degrees(phi) + np.degrees((1 -
            self.e2*(np.sin(phi)**2)**2)/(2*np.cos(phi))*(
            (q/(1-self.e2)) - (np.sin(phi)/(1-self.e2*(np.sin(phi)**2))) +
            (1/(2*self.e))*np.log((1 - self.e*np.sin(phi))/(1 +
            self.e*np.sin(phi))))))

def fastfuels_access(domain_params,out_dir): # NEEDS CLEANING UP # DRG, THIS CAN BE ADDED TO THE DOM CLASS ZACH PRODUCED
    import fastfuels as ff
    x_center,y_center,x_ext,y_ext =  domain_params
    albers = AlbersEqualAreaConic()
    fio = ff.open('https://wifire-data.sdsc.edu:9000/fastfuels/index.fio', ftype='s3')
    fio.cache_limit = 1e16
    
    lat, long = albers.inverse(x_center,y_center)
    roi = fio.query(long, lat,xlen=x_ext,ylen=y_ext)
    roi.write(out_dir+'/',res_xyz=[2,2,1])     
    del roi

class QF_templates:
    def __init__(self,dest_path):
        self.dest_path = dest_path
    
    def generate_all(self):
        self.print_QFire_Advanced_User_Inputs()
        self.print_QFire_Bldg_Advanced_User_Inputs()
        self.print_QFire_Plume_Advanced_User_Inputs()
        self.print_QP_buildout()
        self.print_QUIC_fire()
        self.print_QU_buildings()
        self.print_QU_fileoptions()
        self.print_QU_metparams()
        self.print_QU_movingcoords()
        self.print_QU_simparams()
        self.print_Runtime_Advanced_User_Inputs()
        self.print_sensor1()
        self.print_topo()

    def print_QFire_Advanced_User_Inputs(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QFire_Advanced_User_Inputs.template')
        with open(abs_path,'w') as input_file:
            input_file.write("{{frac_cells_fire|default(0.05)}}			! N/A, fraction of the cells on fire that will launch a firebrand\n")
            input_file.write("{{scl_fac_rad|default(40.0)}}			! N/A, scaling factor of the radius represented by the firebrands launched\n")
            input_file.write("{{ts_firebrands_trj|default(1)}}				! s, time step for the firebrands trajectory calculation\n")
            input_file.write("{{ts_firebrands_launch|default(10)}}				! s, how often to launch firebrands\n")
            input_file.write("{{num_firebrands|default(500)}}			! N/A, number of firebrands distributed over the landing area\n")
            input_file.write("{{fb_fraction_launch|default(20.0)}}			! N/A, FB_FRACTION_LAUNCHED_to_RT_ratio\n")
            input_file.write("{{min_b_value_coef|default(50.0)}}			! N/A, min_b_value_coef\n")
            input_file.write("{{fb_frac_of_max_size|default(0.75)}}			! N/A, fb_frac_of_max_size\n")
            input_file.write("{{germination_delay|default(180)}}				! s, germination_delay\n")
            input_file.write("{{frac_cell_fire|default(10.0)}}				! N/A, fraction of the cell on fire (to scale w)\n")
            input_file.write("{{min_num_ign_firebrands|default(50)}}				! N/A, minimum number of ignitions via firebrands at a point\n")
            input_file.write("{{max_num_ign_firebrands|default(100)}}			! N/A, maximum number of ignitions via firebrands at a point\n")
            input_file.write("{{min_theta_value|default(0.523598)}}		! rad, min_theta_value (pi/6)\n")
       
    def print_QFire_Bldg_Advanced_User_Inputs(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QFire_Bldg_Advanced_User_Inputs.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("{{quic_urb_blds_flag|default(1)}}				! N/A, flag to convert QUIC-URB buildings to fuel (0 = no, 1 = yes)\n")
            input_file.write("{{fuel_density_blds|default(0.5)}}			! kg/m3, thin fuel density within buildings (if no fuel is specified)\n")
            input_file.write("{{att_coef_blds|default(2.0)}}			! N/A, attenuation coefficient within buildings\n")
            input_file.write("{{surf_rgh_blds|default(0.1)}}	  ! m, surface roughness within buildings\n")
            input_file.write("{{conv_fuel_cnp_winds|default(1)}}				! N/A, flag to convert fuel to canopy for winds (0 = no, 1 = yes)\n")
            input_file.write("{{upd_cnp_winds_fuel|default(1)}}				! N/A, update canopy winds when fuel is consumed\n")
            input_file.write("{{att_coef_fuel|default(1.0)}}			! N/A, attenuation coefficient within fuel\n")
            input_file.write("{{srf_rgh_fuel|default(0.1)}}	  ! m, surface roughness within fuel\n")
            input_file.write("\n")

    def print_QFire_Plume_Advanced_User_Inputs(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QFire_Plume_Advanced_User_Inputs.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("{{max_plume_ts|default(50000)}}			! N/A, max number of plume at each time step\n")
            input_file.write("{{min_vert_vel_plume|default(0.5)}}			! m/s, minimum vertical velocity of a plume. If wc is below minimum, the plume is eliminated\n")
            input_file.write("{{max_vert_vel_plume|default(10)}}			! m/s, maximum vertical velocity of a plume\n")
            #RRL conversation on 10/20/2021 for how to calc MSR:
            # if 0.1 * ws > 0.5: ; msr = round(0.5/ws, 3) ; else: msr = 0.1
            input_file.write("{{(0.5/wind_speed)|round(3) if 0.1*wind_speed>0.5 else 0.1}}			! N/A, minimum speed ratio (plume vertical velocity/wind speed). If below, the plume is eliminated\n")
            input_file.write("{{brunt_vaisala_frq_squared|default(0)}}			! 1/s^2, brunt vaisala frequency squared\n")
            input_file.write("{{creeping_flag|default(1)}}			! N/A creeping flag: 0 = off, 1 = on\n")
            input_file.write("{{ts_plume_flag|default(0)}}			! N/A, plume time step flag (0 = fixed, 1 = adaptable)\n")
            input_file.write("{{ts_plume|default(1)}}			! s, time step plume\n")
            input_file.write("{{sor_opt|default(1)}}			! s, sor option\n")
            input_file.write("{{alpha_fire_cells|default(2.0)}}			! N/A, alpha 2 for fire cells\n")
            input_file.write("{{frac_ign_ener_2atmos|default(0)}}			! N/A, fraction of ignition energy that goes into en 2 atmos\n")
            input_file.write("{{max_angle_2plumes|default(30)}}			! deg, max angle for merging two plumes\n")
            input_file.write("{{frac_plume_length|default(0.7)}}			! N/A, fraction of a plume length that a point on a second plume can be for merging\n")
            input_file.write("{{max_height_btms_plumes|default(100.0)}}    !Max height bottoms of plumes before removed (<=0 sets to max domain height)\n")
            input_file.write("\n")
            input_file.write("\n")

    def print_QP_buildout(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QP_buildout.template')
        with open(abs_path,'w') as input_file:
            input_file.write("{{num_blds|default(0)}}  ! total number of buildings\n")
            input_file.write("{{num_veg_cnps|default(0)}}  ! total number of vegitative canopies\n") 

    def print_QUIC_fire(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QUIC_fire.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("{{fire_flag|default(1)}}					! Fire flag: 1 = for fire; 0 = no fire\n")
            input_file.write("{{random|default(-1)}}				! Random number generator: -1: use time and date, any other integer > 0 is used as the seed\n")
            input_file.write("! FIRE TIMES\n")
            input_file.write("{{fire_ign_unix_epoch_time|default(1488794400)}}		! When the fire is ignited in Unix Epoch time (integer seconds since 1970/1/1 00:00:00)\n")
            input_file.write("{{total_time|default(100)}}				! Total simulation time for the fire [s] \n")
            input_file.write("{{time_step|default(1)}}		   		! time step for the fire simulation [s]\n")
            input_file.write("{{ts_wind_field|default(1)}}					! Number of fire time steps done before updating the quic wind field (integer, >= 1)\n")
            input_file.write("{{ts_print_fire_files|default(25)}}					! After how many fire time steps to print out fire-related files (excluding emissions and radiation)\n")
            input_file.write("{{ts_print_wind_files|default(25)}}					! After how many quic updates to print out wind-related files\n")
            input_file.write("{{ts_print_avg_em_rad_files|default(25)}}					! After how many fire time steps to average emissions and radiation\n")
            input_file.write("{{ts_print_avg_wind_files|default(25)}}					! After how many quic updates to print out averaged wind-related files\n")
            input_file.write("! FIRE GRID\n")
            input_file.write("{{num_vertical_layers|default(1)}}					! Number of vertical layers of fire grid cells (integer)\n")
            input_file.write("{{x_fire_grid_ratio|default(1)}}					! x - fire grid ratio = (QUIC-URB cell size)/(fire cell size), integer, can only be >= 1\n")
            input_file.write("{{y_fire_grid_ratio|default(1)}}					! y - fire grid ratio = (QUIC-URB cell size)/(fire cell size), integer, can only be >= 1\n")
            input_file.write("{{vertical_stretching_flag|default(0)}}					! Vertical stretching flag: 0 = uniform dz, 1 = custom\n")
            input_file.write("{% if vertical_stretching_flag==1 %}")
            input_file.write("{% for nz in qf_nz %}")
            input_file.write("{{ nz }}\n")
            input_file.write("{% endfor %}")
            input_file.write("{% else %}")
            input_file.write("{{fire_grid_cell_height|default(1)}}					! Fire grid cell height (m)\n")
            input_file.write("{% endif %}")
            input_file.write("! FOLDER OF TREES AND IGNITION FILES (full path, empty line if none) -- USE FILE SEPARATOR AT THE END\n")
            input_file.write("{{folder_tree_ignition_files|default(\"\")|tojson}}\n")
            input_file.write("{{fuel_loc_flag|default(1)}}			! 1 = all fuels in one file, 2 = separate files\n")
            input_file.write("{{fuel_type_flag|default(1)}}			! 1 = stream, 2 = with headers\n")
            input_file.write("! FUEL\n")
            input_file.write("{{fuel_density_flag|default(5)}}					! fuel density flag: 1 = uniform; 2 = provided thru QF_FuelDensity.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
            input_file.write("{% if fuel_density_flag==1 %}")
            input_file.write("{{fuel_density_val|default(0.7)}}\n")
            input_file.write("{% endif %}")
            input_file.write("{{fuel_moist_flag|default(5)}}					! fuel moisture flag: 1 = uniform; 2 = provided thru QF_FuelMoisture.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
            input_file.write("{% if fuel_moist_flag==1 %}")
            input_file.write("{{fuel_moist_val|default(0.05)}}\n")
            input_file.write("{% endif %}")
            input_file.write("{% if fuel_density_flag==1 %}") # DRG Add another for loop
            input_file.write("{{fuel_height_flag|default(1)}}					! fuel height flag\n")
            input_file.write("{{fuel_height_val|default(0.7)}}\n")
            input_file.write("{% endif %}")
            input_file.write("! IGNITION LOCATIONS\n")
            input_file.write("{{ign_loc_flag|default(7)}}					! 1 = rectangle, 2 = square ring, 3 = circular ring, 4 = file (QF_Ignitions.inp), 5 = time-dependent ignitions (QF_IgnitionPattern.inp), 6 = ignite.dat (firetech)\n")
            input_file.write("{% if ign_loc_flag==1 %}")
            input_file.write("{{ign_xmin|default(50)}}\n")
            input_file.write("{{ign_ymin|default(150)}}\n")
            input_file.write("{{ign_xlen|default(4)}}\n")
            input_file.write("{{ign_ylen|default(100)}}\n")
            input_file.write("{% endif %}")
            input_file.write("{{ign_loc_param1|default(20)}}\n")
            input_file.write("! FIREBRANDS\n")
            input_file.write("{{firebrands_flag|default(0)}}				! 0 = off, 1 = on\n")
            input_file.write("! OUTPUT FILES (formats depend on the grid type flag)\n")
            input_file.write("{{out_grid_energy_atm|default(0)}}					! Output gridded energy-to-atmosphere (fire grid)\n")
            input_file.write("{{out_cmp_reaction_rate|default(0)}}					! Output compressed array reaction rate (fire grid)\n")
            input_file.write("{{out_cmp_fuel_density|default(1)}}					! Output compressed array fuel density (fire grid)\n")
            input_file.write("{{out_grid_wind|default(1)}}					! Output gridded wind (u,v,w,sigma) (fire grid)\n")
            input_file.write("{{out_grid_qu_wind|default(0)}}					! Output gridded QU winds with fire effects, instantaneous (QUIC-URB grid)\n")
            input_file.write("{{out_grid_qu_avg_wind|default(0)}}					! Output gridded QU winds with fire effects, averaged (QUIC-URB grid)\n")
            input_file.write("{{out_plume_traj|default(0)}}					! Output plume trajectories\n")
            input_file.write("{{out_cmp_fuel_moist|default(0)}}					! Output compressed array fuel moisture (fire grid)\n")
            input_file.write("{{out_vert_mass_burnt|default(0)}}					! Output vertically-integrated % mass burnt (fire grid)\n")
            input_file.write("{{out_grid_plume_locs|default(0)}}					! Output gridded file with plumes locations (QUIC-URB grid)\n")
            input_file.write("{{out_cmp_emis|default(0)}}					! Output compressed array emissions (fire grid)\n")
            input_file.write("{{out_grid_thermal_rad|default(0)}}					! Output gridded thermal radiation (fire grid)\n")

    def print_QU_buildings(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QU_buildings.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("!QUIC 6.26\n")
            input_file.write("{{wall_rgh_length|default(0.1)}}			!Wall roughness length (m)\n")
            input_file.write("{{num_blds|default(0)}}	!Number of Buildings\n")
            input_file.write("{{num_poly_blds_nodes|default(0)}}	!Number of Polygon Building Nodes\n")

    def print_QU_fileoptions(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QU_fileoptions.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("!QUIC 6.26\n")
            input_file.write("{{out_file_format_flag|default(2)}}   !output data file format flag (1=ascii, 2=binary, 3=both)\n")
            input_file.write("{{uofield_flag|default(0)}}   !flag to write out non-mass conserved initial field (uofield.dat) (1=write,0=no write)\n")
            input_file.write("{{uosensorfield_flag|default(0)}}   !flag to write out the file uosensorfield.dat, the initial sensor velocity field (1=write,0=no write)\n")
            input_file.write("{{qu_staggered_velocity_flag|default(0)}}   !flag to write out the file QU_staggered_velocity.bin used by QUIC-Pressure(1=write,0=no write)\n")
            input_file.write("{{fire_energy_per_timestep|default(0)}}    !Output fire energy per timestep\n")

    def print_QU_metparams(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QU_metparams.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("!QUIC 6.26\n")
            input_file.write("{{met_input_flag|default(0)}} !Met input flag (0=QUIC,1=WRF,2=ITT MM5,3=HOTMAC)\n")
            input_file.write("{{num_measu_sites|default(1)}} !Number of measuring sites\n")
            input_file.write("{{max_size_data_points|default(1)}} !Maximum size of data points profiles\n")
            input_file.write("{{site_name|default('sensor1')}} !Site Name\n")
            input_file.write("!File name\n")
            input_file.write("{{sensor_file_name|default('sensor1.inp')}}\n")

    def print_QU_movingcoords(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QU_movingcoords.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("!QUIC 6.3\n")
            input_file.write("{{moving_coord_flag|default(0)}}   !Moving coordinates flag (0=no, 1=yes)\n")
            input_file.write("{{ref_bearing_ship_rel|default(0)}}   !Reference bearing of the ship relative to the non-rotated domain (degrees)\n")
            input_file.write("!Time (Unix Time), Ship Speed (m/s), Ship Bearing (deg), Ocean Current Speed (m/s), Ocean Current Direction (deg)\n")
            input_file.write("1488794400	0	0	0	0\n")
            input_file.write("1488794450	0	0	0	0\n")
            input_file.write("1488794500	0	0	0	0\n")
            input_file.write("1488794550	0	0	0	0\n")
            input_file.write("1488794600	0	0	0	0\n")
            input_file.write("1488794650	0	0	0	0\n")
            input_file.write("1488794700	0	0	0	0\n")
            input_file.write("1488794750	0	0	0	0\n")
            input_file.write("1488794800	0	0	0	0\n")
            input_file.write("1488794850	0	0	0	0\n")
            input_file.write("1488794900	0	0	0	0\n")
            input_file.write("1488794950	0	0	0	0\n")
            input_file.write("1488795000	0	0	0	0\n")

    def print_QU_simparams(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'QU_simparams.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("!QUIC 6.26\n")
            input_file.write("{{nx|default(200)}} !nx - Domain Length(X) Grid Cells\n")
            input_file.write("{{ny|default(200)}} !ny - Domain Width(Y) Grid Cells\n")
            input_file.write("{{nz|default(20)}} !nz - Domain Height(Z) Grid Cells\n")
            input_file.write("{{dx|default(2)}} !dx (meters)\n")
            input_file.write("{{dy|default(2)}} !dy (meters)\n")
            input_file.write("{{vert_stretch_flag|default(3)}} !Vertical stretching flag(0=uniform,1=custom,2=parabolic Z,3=parabolic DZ,4=exponential)\n")
            input_file.write("{{surf_dz|default(1.5)}} !Surface dz (meters)\n")
            input_file.write("{{num_unf_surf_cells|default(4)}} !Number of uniform surface cells\n")
            input_file.write("!dz array (meters)\n")
            input_file.write("{% for dz in dz_array %}")
            input_file.write("{{ dz }}\n")
            input_file.write("{% endfor %}")
            input_file.write("{{total_time_inc|default(1)}} !total time increments\n")
            input_file.write("{{utc_conv|default(0)}} !UTC conversion\n")
            input_file.write("!Begining of time step in Unix Epoch time (integer seconds since 1970/1/1 00:00:00)\n")
            input_file.write("{% if len(wind_timesteps)>0 %}")
            input_file.write("{% for wind_time in wind_timesteps %}")
            input_file.write("{{wind_time|default(1488794400)}}\n")
            input_file.write("{% endfor %}")            
            input_file.write("{% else %}")
            input_file.write("{{start_unix_epoch_time|default(1488794400)}}\n")
            input_file.write("{% endif %}")
            input_file.write("{{rooftop_flag|default(2)}} !rooftop flag (0-none, 1-log profile, 2-vortex)\n")
            input_file.write("{{upwind_cavity_flag|default(3)}} !upwind cavity flag (0-none, 1-Rockle, 2-MVP, 3-HMVP)\n")
            input_file.write("{{street_canyon_flag|default(4)}} !street canyon flag (0-none, 1-Roeckle, 2-CPB, 3-exp. param. PKK, 4-Roeckle w/ Fackrel)\n")
            input_file.write("{{street_intsect_flag|default(1)}} !street intersection flag (0-off, 1-on)\n")
            input_file.write("{{wake_flag|default(3)}} !wake flag (0-none, 1-Rockle, 2-Modified Rockle, 3-Area Scaled)\n")
            input_file.write("{{sidewall_flag|default(1)}} !sidewall flag (0-off, 1-on)\n")
            input_file.write("{{canopy_flag|default(2)}} !Canopy flag (1-Cionco w/o wakes, 2-Cionco w/ wakes)\n")
            input_file.write("{{season_flag|default(1)}} !Season flag (1-Summer, 2-Winter, 3-Transition)\n")
            input_file.write("{{max_num_iters|default(10)}} !Maximum number of iterations\n")
            input_file.write("{{omega_relax_val|default(0.9)}} !omega relax value (For sor solver) Must be between 0 and 2\n")
            input_file.write("{{residual_reduc|default(3)}} !Residual Reduction (Orders of Magnitude)\n")
            input_file.write("{{diffusion_alg_flag|default(0)}} !Use Diffusion Algorithm (1 = on)\n")
            input_file.write("{{num_diffusion_iters|default(20)}} !Number of Diffusion iterations\n")
            input_file.write("{{domain_rot_rel_truenorth|default(0)}} !Domain rotation relative to true north (cw = +)\n")
            input_file.write("{{utmx_domain_orign|default(0.0)}}  !UTMX of domain origin (m)\n")
            input_file.write("{{utmy_domain_orign|default(0.0)}}   !UTMY of domain origin (m)\n")
            input_file.write("{{utm_zone|default(1)}} !UTM zone\n")
            input_file.write("{{utm_zone_letter|default(17)}} !UTM zone leter (1=A,2=B,etc.)\n")
            input_file.write("{{quic_cfd_flag|default(0)}} !QUIC-CFD Flag\n")
            input_file.write("{{explo_bld_dmg_flag|default(0)}} !Explosive building damage flag (1 = on)\n")
            input_file.write("{{bld_array_flag|default(0)}} !Building Array Flag (1 = on)\n")

    def print_Runtime_Advanced_User_Inputs(self):
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'Runtime_Advanced_User_Inputs.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("{{param1|default(8)}}\n")

    def print_sensor1(self): # Needs to be adjusted for multiple timesteps
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'sensor1.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("{{site_name|default(sensor1)}} !Site Name\n")
            input_file.write("{{upper_level_flag|default(0)}} !Upper level flag\n")
            input_file.write("{{upper_level_height|default(50)}} !Upper level height\n")
            input_file.write("{{site_coord_flag|default(1)}} !Site Coordinate Flag\n")
            input_file.write("{{x_coord|default(1)}} !X coordinate\n")
            input_file.write("{{y_coord|default(1)}} !Y coordinate\n")
            input_file.write("{{time_unix_epoch|default(1488794400)}} !Time in unix Epoch Time\n")
            input_file.write("{{site_bdry_layer|default(1)}} !site boundary layer\n")
            input_file.write("{{size_zo|default(0.01)}} !site zo\n")
            input_file.write("{{recp_moninobu_length|default(0.0)}} !reciprocal Monin-Obukhov Length(1/m)\n")
            input_file.write("!height, speed, direction\n")
            input_file.write("{{height|default(10)}} {{wind_speed|default(1.0)}} {{wind_dir|default(270.0)}}\n")

    def print_topo(self): # 
        dest_path = self.dest_path
        abs_path = os.path.join(dest_path,'topo.template')
        with open(abs_path, 'w') as input_file:
            input_file.write("!Relative filepath to topo .dat file (ex:: \"../path/to/topo.dat\")\n")
            # input_file.write("\"\"\n")
            input_file.write("""{{topo_relpath|default("topo.dat")}}\n""")
            input_file.write("{{topo_flag|default(0)}}              !Topo flag 0:Flat 1:Gaussian Hill 3:Constant slope with flat section 5:Custom .dat\n")
            input_file.write("{{smooth_flag|default(1)}}              !Smoothing Flag\n")
            input_file.write("{{topo_smth_iters|default(100)}}			   ! Smoothing iterations (Topo)\n")
            input_file.write("{{total_startup_iters|default(1000)}}              !Total startup iterations\n")
            input_file.write("{{startup_restart_iters|default(0)}}            !Startup Iterations Restart Period (set to zero for no restart)\n")
            input_file.write("{{preconditioning|default(3)}}              !Preconditioning\n")

class IgnitionClass_Rectangular: # NEEDS TO BE ADDED TO A FULL IGNITION FUNCTION
    def __init__(self):
        self.xmin = None        # Defined in calc_ignition
        self.ymin = None        # Defined in calc_ignition
        self.xlen = None        # Defined in calc_ignition
        self.ylen = None        # Defined in calc_ignition

    def calc_ignition(self, user_data): # Done
        pwind_len = 200 #  Length of fireline perpendicular to wind direction
        dwind_len = 10  # Length of fireline parallel to wind direction (downwind)
           
        dx = user_data['QU_simparams']['dx']
        dy = user_data['QU_simparams']['dy']
        nx = user_data['QU_simparams']['nx']
        ny = user_data['QU_simparams']['ny']
           
        wind_dir = user_data['sensor1']['wind_dir']

        if ( 225 < wind_dir <= 315):     
            self.xmin = (dx*nx) * 0.2
            self.ymin = (dy*ny - pwind_len) /2.
            self.xlen = dwind_len
            self.ylen = pwind_len
        elif (135 < wind_dir <= 225):   
            self.xmin = (dx*nx-pwind_len) / 2.
            self.ymin = (dy*ny) *0.2
            self.xlen = pwind_len
            self.ylen = dwind_len
        elif (45 < wind_dir <= 135):
            self.xmin = (dx*nx) * 0.8
            self.ymin = (dy*ny - pwind_len) /2.
            self.xlen = dwind_len
            self.ylen = pwind_len
        else:
            self.xmin = (dx*nx-pwind_len) / 2.
            self.ymin = (dy*ny) *0.8
            self.xlen = pwind_len
            self.ylen = dwind_len
    
class QF_base():
    def __init__(self,dest_path,template_path):
        self.dest_path = dest_path
        self.template_path = template_path
        self.file_names = [
            'QFire_Advanced_User_Inputs',
            'QFire_Bldg_Advanced_User_Inputs',
            'QFire_Plume_Advanced_User_Inputs',
            'QP_buildout',
            'QUIC_fire',
            'QU_buildings',
            'QU_fileoptions',
            'QU_metparams',
            'QU_movingcoords',
            'QU_simparams',
            'Runtime_Advanced_User_Inputs',
            'sensor1',
            'topo',
        ]

    def empty_file_dict(self):
        return dict.fromkeys(self.file_names,{})

    def base_gen_input(self,file_name,data):
        dest_path = self.dest_path
        template_file = file_name + '.template'
        inp_file = file_name + '.inp'
        temp_abs_path = os.path.join(self.template_path,template_file)
        inp_abs_path = os.path.join(dest_path,inp_file)
        with open(temp_abs_path) as file:
            template = Template(file.read())
        try:
            template_with_data = template.render(data)
        except TemplateError as err:
            sys.exit(u'Template Error: {}'.format(err.message))

        with open(inp_abs_path, "w") as out_file:
            out_file.write(template_with_data)


        self.print_QU_landuse()

    def generate_all(self, user_data, **kwargs):
        dest_path = self.dest_path
        if kwargs['seed'] is not None:
            user_data['QUIC_fire']['random'] = kwargs['seed']
    
        if kwargs['output'] is not None:
            user_data['QUIC_fire']['ts_print_fire_files'] = kwargs['output']['steps_fire']
            user_data['QUIC_fire']['ts_print_wind_files'] = kwargs['output']['steps_wind']

        #print('ign', kwargs['ignition'])
        if kwargs['ignition'] is not None:
            ign = kwargs['ignition']

            if 'dat' in ign:
                # contents of ignite.dat specified from inputs
                with open(os.path.join(dest_path, 'ignite.dat'), 'w') as f:
                    f.writelines(ign['dat'])
                user_data['QUIC_fire']['ign_loc_flag'] = 7
                user_data['QUIC_fire']['ign_loc_param1'] = ign['perc']
            else:
                raise ValueError('Unknown type of ignition; missing dat contents.')

        else:
            # default rectangle
            Igns = IgnitionClass_Rectangular()
            Igns.calc_ignition(user_data)
            user_data['QUIC_fire']['ign_loc_flag'] = 1
            user_data['QUIC_fire']['ign_xmin'] = Igns.xmin
            user_data['QUIC_fire']['ign_ymin'] = Igns.ymin
            user_data['QUIC_fire']['ign_xlen'] = Igns.xlen
            user_data['QUIC_fire']['ign_ylen'] = Igns.ylen
            user_data['QUIC_fire']['ign_loc_param1'] =  2


        if kwargs['fire_grid'] is not None and \
            kwargs['fire_grid']['vertical_stretching'] is not None:

            user_data['QUIC_fire']['qf_nz'] = kwargs['fire_grid']['vertical_stretching']
            user_data['QUIC_fire']['num_vertical_layers'] = len(user_data['QUIC_fire']['qf_nz'])


        if kwargs['dz'] is not None:
            user_data['QU_simparams']['dz_array'] = kwargs['dz']
            user_data['QU_simparams']['nz'] = len(kwargs['dz'])

        empty_params_dict = self.empty_file_dict()
        params_dict = {**empty_params_dict, **user_data}
        for file, data in params_dict.items():
            self.base_gen_input(file,data)
    
    def print_QU_landuse(self):
        abs_path = os.path.join(self.dest_path,'QU_landuse.inp')
        with open(abs_path,'wb') as input_file:
            # write constant values - QF does not use them. (what to these mean?)
            np.array([1.3563156e-19, 2.9410407e-38, 1.3563156e-19, 7.7097618e-33],
                     dtype=np.float32).tofile(input_file)
            
