import os, sys
import getopt
import xml.etree.ElementTree as ET
from shutil import rmtree

#from __main__ import name

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------

def get_gml_src_proj(path_gml):
    """Reads the source projection system."""
    f = open(path_gml)
    text = f.read()
    #print 'text',text
    text = text.split('xmlns:gml=')[-1]
    gml = text.split('"')[1]
    
    tree = ET.parse(path_gml)
    
    # GET EOP-URL
    root = tree.getroot()
    eop_url = root.tag.split('}')[0].strip('{')
    
    node = tree.find('.//{'+gml+'}Envelope')
    try:
        srs_string = node.get('srsName')
        epsg_no = srs_string.split(':')[-1]
        epsg_code = 'EPSG:{}'.format(epsg_no)
    except:
        print('Warning, gml has no SRS information.')
        epsg_code = None
    
    return epsg_code

# ----------------------------------------------------------------------------------------------------------------------


def get_mask_shp(path_in_gml, path_out, union = False):
    """Transform the gml cloud mask to a shape file."""
    
    if not os.path.isdir(os.path.split(path_out)[0]):
        os.makedirs(os.path.split(path_out)[0])
    
    name = os.path.split(path_in_gml)[-1]

    if name.endswith('gml'):
        #path_out_shp = os.path.join(path_out,name.replace('.gml','.shp'))
                
        src_srs = get_gml_src_proj(path_in_gml)

        trg_srs = 'EPSG:4326'
        
        if src_srs is None:
            return False
            #cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/ogr2ogr -t_srs {} -f "ESRI Shapefile" {} {}'.format(trg_srs, path_out, path_in_gml)
        else:

            if union:
                tmpFP = os.path.join(os.path.split(path_out)[0],'tmp')
                if not os.path.isdir(tmpFP):
                    os.makedirs(tmpFP)
                tmp_out = os.path.join(tmpFP,'tmp.shp')
     
                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/ogr2ogr -s_srs {} -t_srs {} -f "ESRI Shapefile" {} {} '.format(src_srs, trg_srs, tmp_out, path_in_gml)
                os.system(cmd)

                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/ogr2ogr -s_srs {} -t_srs {} -f "ESRI Shapefile" {} {}  -dialect sqlite -sql "SELECT ST_Union(geometry) AS geometry FROM tmp"'.format(trg_srs, trg_srs, path_out, tmp_out)

                os.system(cmd) 
                rmtree(tmpFP)

            else:
                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/ogr2ogr -s_srs {} -t_srs {} -f "ESRI Shapefile" {} {} '.format(src_srs, trg_srs, path_out, path_in_gml)
            
                os.system(cmd)  
                
            return path_out  
    

# ----------------------------------------------------------------------------------------------------------------------

def SentinelCloud(srcFPN, dstFPN):
    path_out_shp = get_mask_shp(srcFPN, dstFPN)

if __name__ == '__main__':
    #path_in, path_out = parameters(sys.argv[1:])   
    #path_out_shp = get_cloud_shp(path_in, path_out)
    srcFPN = '/Volumes/karttur2tb/sentinel/MSI/S2MSI1C/33/WXR/20170612/R008/S2A_MSIL1C_20170612T104021_N0205_R008_T33WXR_20170612T104023.SAFE/GRANULE/L1C_T33WXR_A010301_20170612T104023/QI_DATA/MSK_CLOUDS_B00.gml'
    dstFPN = '/Volumes/karttur2tb/sentinel/MSI/mask/33/WXR/20170612/R008/cloud.shp'
    path_out_shp = get_mask_shp(srcFPN, dstFPN)