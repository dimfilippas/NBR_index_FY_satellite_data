# -*- coding: utf-8 -*-
"""
Created on Wed May  9 12:37:41 2018

1.fy3b_1000 Convert FY3b hd5 L1 1000m files to geotiff
2.fy3b_250  Convert FY3b hd5 L1 0250m files to geotiff

"""

import copy
import glob
import os
import shutil
from geoimread import geoimread
from geoimwrite import geoimwrite
import sys
import numpy as np
from datetime import datetime

def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""

def remove_by_ext(PathOfImages, extension):
    for root, dirnames, filenames in os.walk(PathOfImages):
        for filename in filenames:
            if filename.endswith(extension):
                os.remove(os.path.join(root, filename))

# Convert FY3b hd5 L1 1000m files to geotiff
def fy3b_1000(data_path,ExportPath):
    #for i in Folders:
    # cd to each folder
    i = data_path
    os.chdir(i)
    # find 1000M_L1B.h5 file
    l1b_1000_h5 = ''.join(glob.glob('*MERSI_1000M_L1B.H5'))
    
    # set 1000M_L1B.vrt file
    l1b_1000_vrt = os.path.join(ExportPath, '.'.join(l1b_1000_h5.split('.')[:-1]) + '.vrt')
    # set 1000M_L1B.tif file
    l1b_1000_tif = os.path.join(ExportPath, '.'.join(l1b_1000_h5.split('.')[:-1]) + '.tif')
    # set gdal_translate command and execute
    cmd_tr_1000 = 'gdal_translate -of VRT HDF5:' + os.path.join(i, l1b_1000_h5) + \
                  '://EV_1KM_RefSB ' + l1b_1000_vrt
    res = os.system(cmd_tr_1000)
    
    # check if hd5 has appropriate subsets
    if res == 0:
        # open vrt file
        in_file = open(l1b_1000_vrt, "rt")
        # read the entire file into a string variable
        contents_1000 = in_file.read()
        # find substring that is going to be replaced
        erase = find_between(contents_1000, "</Metadata>", "<VRTRasterBand")
        # set substring
        rplc = """ 
            <Metadata domain="GEOLOCATION">
            <MDI key="LINE_OFFSET">1</MDI>
            <MDI key="LINE_STEP">1</MDI>
            <MDI key="PIXEL_OFFSET">1</MDI>
            <MDI key="PIXEL_STEP">1</MDI>
            <MDI key="SRS">GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,""" + \
               """AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],""" + \
               """UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]</MDI>
               <MDI key="X_BAND">1</MDI>
               <MDI key="X_DATASET">HDF5:""" + os.path.join(i, l1b_1000_h5) + """://Longitude</MDI>
            <MDI key="Y_BAND">1</MDI>
            <MDI key="Y_DATASET">HDF5:""" + os.path.join(i, l1b_1000_h5) + """://Latitude</MDI>
            </Metadata>"""
        # replace substring
        contents_1000 = contents_1000.replace(erase, rplc)
        in_file.close()
        # write string to vrt file
        fvrt = open(l1b_1000_vrt, "w")
        fvrt.write(contents_1000)
        fvrt.close()
        # set gdalwarp command and execute
        cmd_wr_1000 = 'gdalwarp -ot Float32 -geoloc -t_srs EPSG:4326 -overwrite ' + l1b_1000_vrt + ' ' + l1b_1000_tif
        os.system(cmd_wr_1000)
        # Create and move vrt and tif files to results folder
        if os.path.exists(ExportPath) and os.listdir(ExportPath) == []:
            shutil.rmtree(ExportPath)
        elif os.path.exists(ExportPath) and os.listdir(i) != []:
            shutil.move(l1b_1000_tif, os.path.join(ExportPath, l1b_1000_tif.split('/')[-1]))
            shutil.move(l1b_1000_vrt, os.path.join(ExportPath, l1b_1000_vrt.split('/')[-1]))
        else:
            os.mkdir(ExportPath)
            shutil.move(l1b_1000_tif, os.path.join(ExportPath, l1b_1000_tif.split('/')[-1]))
            shutil.move(l1b_1000_vrt, os.path.join(ExportPath, l1b_1000_vrt.split('/')[-1]))
        os.chdir(ExportPath)
        for n in range(1, 16):
            cmd_tr_bands = 'gdal_translate ' + os.path.join(ExportPath, l1b_1000_tif.split('/')[-1]) + ' -b ' + \
                           str(n) + ' ' + '.'.join(
                ((l1b_1000_tif.split('/')[-1]).split('.')[:-1])) + '_EV_1KM_RefSB_b' + str(n + 5) + '.tif'
            os.system(cmd_tr_bands)
        # os.remove(os.path.join(i, 'results', l1b_1000_tif.split('/')[-1]))
        print ExportPath
    else:
        print "Non appropriate file in", i

# Convert FY3b hd5 L1 250m files to geotiff
def fy3b_250(data_path,ExportPath):
    #for j in Folders:
    # cd to each folder
    j = data_path
    os.chdir(j)
    # set names for appropriate bands
    data_list = ['EV_250_Emissive', 'EV_250_RefSB_b1', 'EV_250_RefSB_b2', 'EV_250_RefSB_b3', 'EV_250_RefSB_b4']
    # find 250M_L1B.h5 file        
    l1b_0250_h5 = ''.join(glob.glob('*MERSI_0250M_L1B.H5'))
    # iterate through each band
    for f in data_list:
        # set 1000M_L1B.vrt file
        l1b_0250_vrt = os.path.join(ExportPath, '.'.join(l1b_0250_h5.split('.')[:-1]) + '_' + f + '.vrt')
        # set 1000M_L1B.tif file
        l1b_0250_tif = os.path.join(ExportPath, '.'.join(l1b_0250_h5.split('.')[:-1]) + '_' + f + '.tif')
        # set gdal_translate command and execute
        cmd_tr_250 = 'gdal_translate -of VRT HDF5:' + os.path.join(j, l1b_0250_h5) + \
                     '://' + f + ' ' + l1b_0250_vrt
        res = os.system(cmd_tr_250)
        print res
        # check if hd5 has appropriate subsets
        if res == 0:
            # open vrt file
            in_file = open(l1b_0250_vrt, "rt")
            # read the entire file into a string variable
            contents_0250 = in_file.read()
            # find substring that is going to be replaced
            erase = find_between(contents_0250, "</Metadata>", "<VRTRasterBand")
            # set substring
            rplc = """ 
                    <Metadata domain="GEOLOCATION">
                    <MDI key="LINE_OFFSET">1</MDI>
                    <MDI key="LINE_STEP">4</MDI>
                    <MDI key="PIXEL_OFFSET">1</MDI>
                    <MDI key="PIXEL_STEP">4</MDI>
                    <MDI key="SRS">GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,""" + \
                   """AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,""" + \
                   """AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]</MDI>
                   <MDI key="X_BAND">1</MDI>
                   <MDI key="X_DATASET">HDF5:""" + os.path.join(j, l1b_0250_h5) + """://Longitude</MDI>
                    <MDI key="Y_BAND">1</MDI>
                    <MDI key="Y_DATASET">HDF5:""" + os.path.join(j, l1b_0250_h5) + """://Latitude</MDI>
                    </Metadata>"""
            # replace substring
            contents_0250 = contents_0250.replace(erase, rplc)
            in_file.close()
            # write string to vrt file
            fvrt = open(l1b_0250_vrt, "w")
            fvrt.write(contents_0250)
            fvrt.close()
            # set gdalwarp command and execute
            cmd_wr_0250 = 'gdalwarp -ot Float32 -geoloc -t_srs EPSG:4326 -overwrite ' + l1b_0250_vrt + ' ' + l1b_0250_tif
            os.system(cmd_wr_0250)
            # Create and move vrt and tif files to results folder
            if os.path.exists(ExportPath) and os.listdir(ExportPath) == []:
                shutil.rmtree(ExportPath)
            elif os.path.exists(ExportPath) and os.listdir(ExportPath) != []:
                shutil.move(l1b_0250_tif, os.path.join(ExportPath, l1b_0250_tif.split('/')[-1]))
                shutil.move(l1b_0250_vrt, os.path.join(ExportPath, l1b_0250_vrt.split('/')[-1]))
            else:
                os.mkdir(ExportPath)
                shutil.move(l1b_0250_tif, os.path.join(ExportPath, l1b_0250_tif.split('/')[-1]))
                shutil.move(l1b_0250_vrt, os.path.join(ExportPath, l1b_0250_vrt.split('/')[-1]))
            res_prn = ExportPath
        else:
            res_prn = "Non appropriate file in" + j
    print res_prn

def reproject(master_path,master_reproject):
    
    os.chdir(master_path)
    file7 = ''.join(glob.glob('*_b7.tif'))
    cmd_file7 = 'gdalwarp -t_srs  "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" ' + file7 + ' ' + master_reproject + '/' + file7
    os.system(cmd_file7)
    
    file4 = ''.join(glob.glob('*_b4.tif'))
    cmd_file4 = 'gdalwarp -t_srs  "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" ' + file4 + ' ' + master_reproject + '/' + file4
    os.system(cmd_file4)
    '''
    fileCMD = ''.join(glob.glob('*_CLM.tif'))
    cmd_fileCMD = 'gdalwarp -t_srs  "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" ' + fileCMD + ' ' + master_reproject + '/' + fileCMD
    os.system(cmd_fileCMD)
    '''

def nbr_index(k,mr,FYobj):
    os.chdir(k)
    #channel_16_name = ''.join(glob.glob('*MERSI_1000M_L1B_EV_1KM_RefSB_b16.tif'))
    #channel_07_name = ''.join(glob.glob('*MERSI_1000M_L1B_EV_1KM_RefSB_b7.tif'))
    FYobj.channel_04_name = ''.join(glob.glob('*_b4.tif'))
    FYobj.channel_04_name_clip = FYobj.channel_04_name[:-4]+'_clip.tif'
    FYobj.channel_07_name = ''.join(glob.glob('*MERSI_1000M_L1B_EV_1KM_RefSB_b7.tif'))
    FYobj.channel_07_resm = FYobj.channel_07_name[:-4] + '_resampled.tif'
    FYobj.channel_07_r_cl = FYobj.channel_07_name[:-4] + '_resampled_clip.tif'
    FYobj.channel_CLM     =''.join(glob.glob('*_CLM*.tif'))
    
    cmd_gdalwrap = 'gdalwarp -r lanczos %s %s %s' %(FYobj.channel_04_name,FYobj.channel_07_name,os.path.join(k,FYobj.channel_07_resm))
    os.system(cmd_gdalwrap)
    
    
    cmd_clip_04   = 'gdalwarp -srcnodata 0 -ot Float32 -q -of GTiff -cutline %s -crop_to_cutline -of GTiff %s %s -overwrite' %(FYobj.aoi,FYobj.channel_04_name,FYobj.channel_04_name_clip)
    os.system(cmd_clip_04)
    cmd_clip_07   = 'gdalwarp -srcnodata 0 -ot Float32 -q -of GTiff -cutline %s -crop_to_cutline -of GTiff %s %s -overwrite' %(FYobj.aoi,FYobj.channel_07_resm,FYobj.channel_07_r_cl)
    os.system(cmd_clip_07)
    channel_04  = geoimread(os.path.join(k, FYobj.channel_04_name_clip))
    channel_07  = geoimread(os.path.join(k, FYobj.channel_07_r_cl))
    
    nbr = (channel_04[0]-channel_07[0])/(channel_04[0] + channel_07[0])
    
    geoimwrite(os.path.join('NBR' + mr + '.tif'), nbr, channel_04[1], channel_04[2], channel_04[3])

def dnbr_index(FYobj):
    os.chdir(FYobj.master_reproject)
    FYobj.nbr_post_path = ''.join(glob.glob('NBRmaster.tif'))
    nbr_post      = geoimread(FYobj.nbr_post_path)
    os.chdir(FYobj.reference_reproject)
    FYobj.nbr_pre_path  = ''.join(glob.glob('NBRreference.tif'))
    
    cmd_align = 'gdalwarp %s %s %s' %(os.path.join(FYobj.master_reproject,FYobj.nbr_post_path),os.path.join(FYobj.reference_reproject,FYobj.nbr_pre_path),os.path.join(FYobj.master_reproject,'NBRreference_aligned.tif'))
    os.system(cmd_align)
    nbr_pre_aligned = geoimread(os.path.join(FYobj.master_reproject,'NBRreference_aligned.tif'))
    cmd_clip_master    = 'gdalwarp -ot Float32 -q -of GTiff -cutline %s -crop_to_cutline -of GTiff %s %s -overwrite' %(FYobj.aoi,os.path.join(FYobj.master_reproject,FYobj.nbr_post_path),os.path.join(FYobj.results_path,'NBRmaster.tif'))
    cmd_clip_reference = 'gdalwarp -ot Float32 -q -of GTiff -cutline %s -crop_to_cutline -of GTiff %s %s -overwrite' %(FYobj.aoi,os.path.join(FYobj.master_reproject,'NBRreference_aligned.tif'),os.path.join(FYobj.results_path,'NBRreference.tif'))
    os.system(cmd_clip_master)
    os.system(cmd_clip_reference)
    
    os.chdir(FYobj.master_reproject)
    FYobj.CLM_master = ''.join(glob.glob('*_CLM_clip.tif'))
    channel_CLM_m = geoimread(os.path.join(FYobj.master_reproject, FYobj.CLM_master))
    
    os.chdir(FYobj.reference_reproject)
    FYobj.CLM_reference = ''.join(glob.glob('*_CLM_clip.tif'))
    
    cmd_wrap_cld = 'gdalwarp -q -of GTiff -cutline %s -crop_to_cutline -overwrite %s %s %s' %(FYobj.aoi,os.path.join(FYobj.master_reproject,'NBRreference_aligned.tif'),FYobj.CLM_reference,FYobj.CLM_reference)
    os.system(cmd_wrap_cld)
    
    channel_CLM_r = geoimread(os.path.join(FYobj.reference_reproject, FYobj.CLM_reference))
    
    os.chdir(FYobj.results_path)
    nbr_master    = geoimread(os.path.join(FYobj.results_path, 'NBRmaster.tif'))
    nbr_reference = geoimread(os.path.join(FYobj.results_path, 'NBRreference.tif'))
    
    dnbr = nbr_reference[0]-nbr_master[0]
    
    try:
        dnbr[channel_CLM_m[0]>0]=-9999
    except:
        pass
    
    try:
        dnbr[channel_CLM_r[0]>0]=-9999
    except:
        pass
    
    geoimwrite(FYobj.results_path + '/dNBR.tif', dnbr, nbr_master[1], nbr_master[2], nbr_master[3])
    
    cmd_nodata = 'gdalwarp -srcnodata -9999 dNBR.tif dNBR_nodata.tif'
    os.system(cmd_nodata)
    
def CloudMask_FY3b(image_path,FYobj,mr):
    # Create Cloud Mask for FY-3b
    # Cloud Mask has three classes:
    #        0 -> Non Cloud
    #        1 -> Propably Cloud
    #        2 -> Cloud
    
    os.chdir(image_path)
    b1_tif = ''.join(glob.glob('*L1B_EV_250_RefSB_b1.tif'))
    b2_tif = ''.join(glob.glob('*L1B_EV_250_RefSB_b2.tif'))
    b3_tif = ''.join(glob.glob('*L1B_EV_250_RefSB_b3.tif'))
    b1 = geoimread( b1_tif)
    b2 = geoimread( b2_tif)
    b3 = geoimread( b3_tif)

    d = np.dstack([b1[0], b2[0], b3[0]])

    low_values_flags = d > 900

    d[low_values_flags] = 1

    d[~low_values_flags] = 0

    clm = copy.copy(b1[0])

    clm[:, :] = 0

    mask_cloud = np.all(d == (1, 1, 1), axis=-1)

    mask_nc_1 = np.all(d == (0, 0, 0), axis=-1)

    mask_nc_234 = np.logical_or(np.all(d == (0, 1, 0), axis=-1),
                                np.all(d == (0, 0, 1), axis=-1),
                                np.all(d == (0, 1, 1), axis=-1))

    mask_non_cloud = np.logical_or(mask_nc_1, mask_nc_234)

    mask_pc = np.logical_or(np.all(d == (1, 0, 0), axis=-1),
                            np.all(d == (0, 1, 0), axis=-1),
                            np.all(d == (1, 0, 1), axis=-1))
    clm[mask_non_cloud] = 0
    clm[mask_pc]        = 1
    clm[mask_cloud]     = 2
    geoimwrite(b1_tif[:-7]+'_CLM.tif', clm, b1[1], b1[2], b1[3])
    if mr == 'master':
        cmd_clip_clm   = 'gdalwarp -srcnodata 0 -t_srs  "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" -ot Float32 -q -of GTiff -cutline %s -crop_to_cutline -of GTiff %s %s -overwrite' %(FYobj.aoi,b1_tif[:-7]+'_CLM.tif',os.path.join(FYobj.master_reproject,b1_tif[:-7]+'_CLM_clip.tif'))
    else:
        cmd_clip_clm   = 'gdalwarp -srcnodata 0 -t_srs  "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" -ot Float32 -q -of GTiff -cutline %s -crop_to_cutline -of GTiff %s %s -overwrite' %(FYobj.aoi,b1_tif[:-7]+'_CLM.tif',os.path.join(FYobj.reference_reproject,b1_tif[:-7]+'_CLM_clip.tif'))        
    os.system(cmd_clip_clm)

def CreatePaths(FYobj):
    if not os.path.exists(FYobj.ExportPath):
        os.makedirs(FYobj.ExportPath)
        
    if not os.path.exists(FYobj.master_path):
        os.makedirs(FYobj.master_path)
        
    if not os.path.exists(FYobj.reference_path):
        os.makedirs(FYobj.reference_path)    
    
    if not os.path.exists(FYobj.master_reproject):
        os.makedirs(FYobj.master_reproject)
        
    if not os.path.exists(FYobj.reference_reproject):
        os.makedirs(FYobj.reference_reproject)
        
    if not os.path.exists(FYobj.results_path):
        os.makedirs(FYobj.results_path)

def Clean(FYobj):
    for CleanUp in glob.glob(FYobj.results_path + '/*.*'):
        if not CleanUp.endswith('dNBR_nodata.tif'):    
            os.remove(CleanUp)
    os.rename('dNBR_nodata.tif', 'dNBR.tif')
    
def main():
    startTime = datetime.now()
    # Set Path of FY3b images
    FYobj = lambda: None
    FYobj.ExportPath          = '/home/df/Projects/NBR_FY3/nbr'
    FYobj.master_image        = '/mnt/Xband_Storage/X/FY3/20180806_1335/'
    FYobj.reference_image     = '/mnt/Xband_Storage/X/FY3/20180718_1321/'
    FYobj.aoi                 = '/home/df/Projects/NBR_FY3/Aux_data/greeceaoiv_sinu.shp'
    FYobj.master_path         = os.path.join(FYobj.ExportPath,'master')
    FYobj.reference_path      = os.path.join(FYobj.ExportPath,'reference')
    FYobj.master_reproject    = os.path.join(FYobj.ExportPath,'master_reprojected')
    FYobj.reference_reproject = os.path.join(FYobj.ExportPath,'reference_reprojected')
    FYobj.results_path        = os.path.join(FYobj.ExportPath,'results')
    FYobj.master              = 'master'
    FYobj.reference           = 'reference'
    
    CreatePaths(FYobj)
    
    # Create geotiff
    
    fy3b_1000(FYobj.master_image,FYobj.master_path)
    fy3b_250(FYobj.master_image,FYobj.master_path)
    fy3b_1000(FYobj.reference_image,FYobj.reference_path)
    fy3b_250(FYobj.reference_image,FYobj.reference_path)
    '''
    CloudMask_FY3b(FYobj.master_path,FYobj,FYobj.master)
    CloudMask_FY3b(FYobj.reference_path,FYobj,FYobj.reference)
    
    reproject(FYobj.master_path,FYobj.master_reproject)
    reproject(FYobj.reference_path,FYobj.reference_reproject)
    
    nbr_index(FYobj.master_reproject,FYobj.master,FYobj)
    nbr_index(FYobj.reference_reproject,FYobj.reference,FYobj)
    
    dnbr_index(FYobj)
    Clean(FYobj)
    # Remove xml files
    extension = ['.xml', '.vrt']
    remove_by_ext(FYobj.ExportPath, extension[0])
    remove_by_ext(FYobj.ExportPath, extension[1])
    print datetime.now() - startTime
    '''
if __name__ == '__main__':
    sys.exit(main())
