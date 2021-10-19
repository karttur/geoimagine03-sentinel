'''
Created on 8 juni 2018

@author: thomasgumbricht
'''

import geoimagine.gis.mj_gis_v80 as mj_gis
import geoimagine.zipper.explode as zipper
import geoimagine.support.karttur_dt as mj_dt
#from geoimagine.kartturmain import RasterProcess
#import geoimagine.sentinel.gml_transform as gml_transform
from geoimagine.kartturmain import Composition, LayerCommon, VectorLayer, RasterLayer
#from sentinelsat.sentinel import SentinelAPI
from sentinelsat import SentinelAPI
import os
import xml.etree.ElementTree as ET
from shutil import rmtree, move
import subprocess

class SentinelComposition:
    '''
    class for sentinel compositions
    '''
    def __init__(self, compD):  
        for key in compD:
            if '_' in compD[key]:
                exitstr = 'the "%s" parameter can not contain underscore (_): %s ' %(key, compD[key])
                exit(exitstr) 
            setattr(self, key, compD[key])
        if not hasattr(self, 'folder'):
            exitstr = 'All sentinel compositions must contain a folder'
            exit(exitstr)
    
class SentinelTile(LayerCommon):
    '''Class for sentinel tiles'''
    def __init__(self, composition, locusD, datumD, filepath, FN): 
        """The constructor expects an instance of the composition class."""
        LayerCommon.__init__(self)

        self.comp = composition
        
        self.locus = locusD['locus']
        self.locuspath = locusD['path']

        self.path = filepath
        self.FN = FN

        self.datum = lambda: None
        for key, value in datumD.items():
            setattr(self.datum, key, value)
        if self.datum.acqdate:
            self._SetDOY()
            self._SetAcqdateDOY()
        self._SetPath()
        
    def _SetPath(self):
        """Sets the complete path to sentinel tiles"""
        
        self.FP = os.path.join('/Volumes',self.path.volume, self.comp.system, self.comp.source, self.comp.division, self.comp.folder, self.locuspath, self.datum.acqdatestr)
        self.FPN = os.path.join(self.FP,self.FN)
        if ' ' in self.FPN:
            exitstr = 'EXITING sentinel FPN contains space %s' %(self.FPN)
            exit(exitstr)
     
class ProcessSentinel:
    '''class for senitnel specific processing'''   
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        print (self.process.proc.processid) 

        #direct to subprocess
        if self.process.proc.processid == 'downloadSentinelTile':
            self._DownloadSentinelTile(self.process.params.tileid)
        elif 'downloadSentinelData' in self.process.proc.processid:
            self._DownloadSentinelData()
        elif self.process.proc.processid[0:8] == 'download':
            self._LoopSearchLayers()
        elif self.process.proc.processid == 'explodeSentinel':
            self._ExplodeSentinel()
        elif self.process.proc.processid == 'extractsentineltilecoords':
            self._ExtractSentinelTileCoords('sentinel')
        elif self.process.proc.processid == 'extractmgrscoords':
            self._ExtractSentinelTileCoords('mgrs')
        elif self.process.proc.processid == 'LinkDefaultRegionsToSentinel':
            self._LinkDefaultRegionsToSentinel()
        elif self.process.proc.processid == 'findgranuletiles':
            self._FindGranuleTiles()
        elif self.process.proc.processid == 'reorganisesentinel':
            self._ReOrganiseSentinel()    
        elif self.process.proc.processid == 'geochecksentineltiles':
            self._GeoCheckSentinelTiles()
        else:
            exitStr = 'EXITING, unknown process in ProcessSentine (sentinel.sentinel.py: %(s)s)' %{'s':self.process.proc.processid}
            exit(exitStr)
            
    def _ExtractSentinelTileCoords(self, mgrsversion):
        '''Extracts MGRS tile corners from MGRS polygons
        '''
        for locus in self.process.srcLayerD:
            if len(self.process.srcLayerD[locus]) == 0:
                exitstr = 'EXITING, no dates defined in Sentinel._ExtractSentinelTileCoords'
                exit(exitstr)
            for datum in self.process.srcLayerD[locus]:
                if len(self.process.srcLayerD[locus][datum]) == 0:
                    exitstr = 'EXITING, no compositions defined in Sentinel._ExtractSentinelTileCoords'
                    exit(exitstr)
                for comp in self.process.srcLayerD[locus][datum]:
                    self.srcLayer = self.process.srcLayerD[locus][datum][comp]
                    if mgrsversion == 'sentinel':
                        self._ExtractVectorCoords(self.srcLayer.FPN)
                    else:
                        self._ExtractMGRSCoords(self.srcLayer.FPN)
                    
    def _ExtractVectorCoords(self,srcFPN):
        '''Open layer and loop over features
        '''
        southL = ['C','D','E','F','G','H','J','K','L','M']
        northL = ['N','P','Q','R','S','T','U','V','W','X']
        from operator import itemgetter
        #from math import sqrt
        srcDS,srcLayer = mj_gis.ESRIOpenGetLayer(srcFPN)

        fieldName = srcLayer.fieldDefL[0].name
        #Loop over the features in the layer
        for feature in srcLayer.layer:
            mgrs = feature.GetField(fieldName)
            #utmzone = int(mgrs[0:2])
            latzone = mgrs[2:3]
            printstr = 'Extracting Sentinel MGRS %(m)s' %{'m':mgrs}
            print (printstr)
            if mgrs[0:2] in ['01','59','60']:
                #Zone 01 wraps the dateline, and I did not solve that yet, but very few tiles are affected for my stuff at the moment
                continue
            #Craete a geom instance
            geom = mj_gis.Geometry()
            #add the feature and extract the geom
            geom.GeomFromFeature(feature)
            if srcLayer.geomtype.lower() != 'polygon':
                exit('_ExtractVectorCoords can nonly have polygon as input')
                #Insert the centroid of each mgrs zone
                query = {'mgrs':mgrs,'utmzone':int(mgrs[0:2]), 'mgrsid':mgrs[2:6],'centerlon':geom.shapelyGeom.x, 'centerlat':geom.shapelyGeom.y}   
                self.session._InsertTileCoords(query)
            else:
                #Set west, south, east, north from the bounding box
                west, south, east, north = geom.shapelyGeom.bounds

                #Get the coords of the original polygon
                coords = list(geom.shapelyGeom.exterior.coords)
                #Pop the fifth coordinate (identicial to the first)
                coords.pop(4)
                
                #Sort the coordinates using y (latitude)
                ysortedCoords = sorted(coords, key=itemgetter(1))
                #After sorting the two first coordinates are the south edge, and the following two the north edge
                minyCoords = ysortedCoords[0:2]
                maxyCoords = ysortedCoords[2:4]
                #From the sorted southern and northern edge coordinates, extract the corners
                cornerD = {}
                if minyCoords[0][0] < minyCoords[1][0]:
                    cornerD['ll'] = minyCoords[0]
                    cornerD['lr'] = minyCoords[1]
                else:
                    cornerD['ll'] = minyCoords[1]
                    cornerD['lr'] = minyCoords[0]               
                if maxyCoords[0][0] < maxyCoords[1][0]:
                    cornerD['ul'] = maxyCoords[0]
                    cornerD['ur'] = maxyCoords[1]
                else:
                    cornerD['ul'] = maxyCoords[1]
                    cornerD['ur'] = maxyCoords[0]
                #Compose a query (dict) 
                query = {'mgrs':mgrs,'utmzone':int(mgrs[0:2]), 'mgrsid':mgrs[2:6],
                         'east':east, 'west':west, 'north':north, 'south':south,
                         'centerlon':geom.shapelyGeom.centroid.x,'centerlat':geom.shapelyGeom.centroid.y}
                query['ullat'] = cornerD['ul'][1]
                query['ullon'] = cornerD['ul'][0]
                query['urlat'] = cornerD['ur'][1]
                query['urlon'] = cornerD['ur'][0]
                query['lrlat'] = cornerD['lr'][1]
                query['lrlon'] = cornerD['lr'][0]
                query['lllat'] = cornerD['ll'][1]
                query['lllon'] = cornerD['ll'][0]
                llPtL = ( cornerD['ul'],cornerD['ur'],cornerD['lr'],cornerD['ll'] )
                lonlatgeom = mj_gis.Geometry()
                lonlatgeom.PointsToMultiPointGeom(llPtL)
                
                # Set the EPSG code for this utm zone
                if query['utmzone'] < 10:
                    utm_band = '0%(utm)d' %{'utm':query['utmzone']}
                else:
                    utm_band = '%(utm)d' %{'utm':query['utmzone']}
                if latzone in northL:
                    epsg_code = '326' + utm_band
                elif latzone in southL:
                    epsg_code = '327' + utm_band
                else:
                    print (query)
                    print (mgrs)
                    exit('Error in _ExtractVectorCoords')
                    
                dstepsg = int(epsg_code)
                query['EPSG'] = dstepsg

                #Set the projections, source projection is always lonlat (4326) dst the UTM zone
                srcproj = mj_gis.MjProj()
                srcproj.SetFromEPSG(4326)
                dstproj = mj_gis.MjProj()
                dstproj.SetFromEPSG(dstepsg)
                utmgeom =srcproj.ReprojectGeom(lonlatgeom, dstproj)
                
                cornerL = [(p.x, p.y) for p in utmgeom.shapelyGeom]
                
                if int(round(cornerL[0][0])) != int(round(cornerL[3][0])):
                    exitstr = 'ERROR in minx', int(round(cornerL[0][0])), int(round(cornerL[3][0]))
                    print ('cornerL',cornerL)
                    print ('mgrs',mgrs)
                    print ('utmzone',query['utmzone'])
                    exit(exitstr)
                else:
                    minx = int(round(cornerL[0][0]))
                    
                if int(round(cornerL[1][0])) != int(round(cornerL[2][0])):
                    exitstr = 'ERROR in maxx', int(round(cornerL[1][0])), int(round(cornerL[2][0]))
                    print ('cornerL',cornerL)
                    print ('mgrs',mgrs)
                    print ('utmzone',query['utmzone'])
                    exit (exitstr)
                else:
                    maxx = int(round(cornerL[1][0]))
                    
                if int(round(cornerL[2][1])) != int(round(cornerL[3][1])):
                    exitstr = 'ERROR in miny', int(round(cornerL[2][1])), int(round(cornerL[3][1]))
                    print ('cornerL',cornerL)
                    print ('mgrs',mgrs)
                    print ('utmzone',query['utmzone'])
                    exit (exitstr)
                else:
                    miny = int(round(cornerL[2][1]))
                    
                if int(round(cornerL[0][1])) != int(round(cornerL[1][1])):
                    exitstr = 'ERROR in maxy', int(round(cornerL[0][1])), int(round(cornerL[1][1]))
                    print ('cornerL',cornerL)
                    print ('mgrs',mgrs)
                    print ('utmzone',query['utmzone'])
                    print ('llPtL',llPtL)
                    print ('coords',coords)
                    exit (exitstr)
                else:
                    maxy = int(round(cornerL[0][1]))

                query['minx'] = minx
                query['miny'] = miny
                query['maxy'] = maxy
                query['maxx'] = maxx
                
                refsize = 20
                refcols = int((maxx - minx) /refsize)
                reflins = int((maxy - miny)/refsize)
                if refcols != reflins or reflins != 5490:
                    exitstr = 'ERROR in lins/cols',reflins,refcols
                    print ('cornerL',cornerL)
                    print ('mgrs',mgrs)
                    print ('utmzone',query['utmzone'])
                    print ('llPtL',llPtL)
                    print ('coords',coords)
                    exit (exitstr)
                query['refsize'] = refsize
                query['refcols'] = refcols
                query['reflins'] = reflins
                
                #Get the proj4 code
                dstproj.SetProj4()
                query['proj4'] = dstproj.proj4
                #Register the tile and its coordinates in the database
                self.session._InsertTileCoords(query)
  
    def _ExtractMGRSCoords(self,srcFPN):
        '''Open layer and loop over features
        '''
        #from operator import itemgetter
        #from math import sqrt
        southL = ['C','D','E','F','G','H','J','K','L','M']
        northL = ['N','P','Q','R','S','T','U','V','W','X']
   
        srcDS,srcLayer = mj_gis.ESRIOpenGetLayer(srcFPN)
        
        fieldName = srcLayer.fieldDefL[0].name
        #Loop over the features in the layer
        for feature in srcLayer.layer:
            mgrs = feature.GetField(fieldName)
            print (mgrs)
            #Check if this mgrs is included in the sentinel tiling system
            query = {'mgrs':mgrs}
            sentinelFlag = self.session._SelectSentinelTile(query)
            printstr = 'Extracting MGRS %(m)s' %{'m':mgrs}
            print (printstr)
            latzone = mgrs[2:3]
            if mgrs[0:2] == '01':
                #Zone 01 wraps the dateline, and I did not solve that yet, but very few tiles are affected for my stuff at the moment
                continue
            #Craete a geom instance
            geom = mj_gis.Geometry()
            #add the feature and extract the geom
            geom.GeomFromFeature(feature)

            #Set west, south, east, north from the bounding box
            west, south, east, north = geom.shapelyGeom.bounds
            #Get the coords of the original polygon
            #coords = list(geom.shapelyGeom.exterior.coords)
            query = {'mgrs':mgrs,'east':east, 'west':west, 'north':north, 'south':south}
            utmzone = int(mgrs[0:2])
            # Set the EPSG code for this utm zone

            if utmzone < 10:
                utm_band = '0%(utm)d' %{'utm':utmzone}
            else:
                utm_band = '%(utm)d' %{'utm':utmzone}
                            
            if latzone in northL:
                epsg_code = '326' + utm_band
            elif latzone in southL:
                epsg_code = '327' + utm_band
            else:
                print (query)
                print (mgrs)
                exit('Error in _ExtractVectorCoords')
                    
                    
            if south < 0:
                epsg_code = '326' + utm_band
            elif north > 0:
                epsg_code = '327' + utm_band
            else:
                print ('north, south',north, south)
                exit('Error in _ExtractMGRSCoords')
            dstepsg = int(epsg_code)

            #Set the projections, source is always lonlat (4326) dst the UTM zone
            srcproj = mj_gis.MjProj()
            srcproj.SetFromEPSG(4326)
            dstproj = mj_gis.MjProj()
            dstproj.SetFromEPSG(dstepsg)
            utmgeom =srcproj.ReprojectGeom(geom, dstproj)
            
            #Simplified squared MGRS geom
            utmxyL = list(utmgeom.shapelyGeom.exterior.coords)

            ptX = [p[0] for p in utmxyL]
            ptY = [p[1] for p in utmxyL]
            minx = int(round(min(ptX) ))
            maxx = int(round(max(ptX) ))
            miny = int(round(min(ptY) ))
            maxy = int(round(max(ptY) ))
            
            query['minx'] = minx
            query['miny'] = miny
            query['maxy'] = maxy
            query['maxx'] = maxx
            if sentinelFlag:
                query['sentinel'] = 'Y'
            else:
                query['sentinel'] = 'N'
            self.session._InsertMGRSCoords(query)
                
    def _LoopSearchLayers(self):
        '''Loop over layers with polygons to search 
        '''
        #Get the fields to save to the meta database
        self.metaTranslator = self.session._GetMetaTranslator()
        #loop over all input files
        for locus in self.process.srcLayerD:
            if len(self.process.srcLayerD[locus]) == 0:
                exitstr = 'EXITING, no dates defined in Sentinel._downloadSentinelFromVector'
                exit(exitstr)
            for datum in self.process.srcLayerD[locus]:
                if len(self.process.srcLayerD[locus][datum]) == 0:
                    exitstr = 'EXITING, no compositions defined in Sentinel._downloadSentinelFromVector'
                    exit(exitstr)
                for comp in self.process.srcLayerD[locus][datum]:
                    self.srcLayer = self.process.srcLayerD[locus][datum][comp]
                    self._SearchFromVector(self.srcLayer.FPN)
                   
    def _SearchFromVector(self,srcFPN):
        '''Open vector layer and loop over seearch features
        '''
        #print (self.process.params.tablesearchid)
        searchid = self.process.params.tablesearchid
        #Open the datasource and get the layer
        srcDS,srcLayer = mj_gis.ESRIOpenGetLayer(srcFPN)
        srcDS.CloseDS()
        self.srcFPN = srcFPN
        #Loop over the features in the layer
        for feature in srcLayer.layer:
            self.searchid = feature.GetField(searchid)

            #Craete a geom instance
            geom = mj_gis.Geometry()
            #add the feature and extract the geom
            geom.GeomFromFeature(feature)
            boundaryPoly = geom.BoundsToPoly()
            #Inititate the copernicus search
            self._CopernicusSearch(boundaryPoly)
       
    def _CopernicusSearch(self,wktGeom):
        '''Search online sentinel tiles
        '''        
        if not mj_dt.IsDatetime(self.process.params.startdate):
            startdate = mj_dt.yyyymmddDate(self.process.params.startdate)
            enddate = mj_dt.yyyymmddDate(self.process.params.enddate)
            
        searchD = self._CheckPreviousSearches(startdate,enddate,self.process.params.platformname,self.process.params.prodtype,self.process.params.cloudmax)

        if not searchD: 
            return

        api = SentinelAPI('thomas.gumbricht', 'iHg-98G-gug-34t', 'https://scihub.copernicus.eu/dhus')
        if searchD['platformname'] == 'Sentinel-1':
            products = api.query(wktGeom,
                                 date = (searchD['startdate'], searchD['enddate']),
                                 platformname = searchD['platformname'],
                                 producttype = searchD['product']) 
            
        elif searchD['platformname'] == 'Sentinel-2':
            products = api.query(wktGeom,
                                 date = (searchD['startdate'], searchD['enddate']),
                                 platformname = searchD['platformname'],
                                 producttype = searchD['product'],
                                 cloudcoverpercentage = (0,searchD['cloudcover']))  
        else:
            exitstr = 'unrecognized platform, %(p)s, in _CopernicusSearch'
            exit(exitstr)
        #Construct empty dictionaries for holding the metadata 
        self.metaquery = {}
        self.tilequery = {}
        #Loop over the metaadata as retrieved from the product search
        for key in products:
            for pos in products[key]:
                #print ('key',pos,products[key][pos])
                if pos in self.metaTranslator and self.metaTranslator[pos]['tab'] in ['meta','metatiles']:
                    if self.metaTranslator[pos]['typ'] == 'integer':
                        self.metaquery[self.metaTranslator[pos]['dst']] = int(products[key][pos])
                    elif self.metaTranslator[pos]['typ'] == 'real':
                        self.metaquery[self.metaTranslator[pos]['dst']] = float(products[key][pos])
                    elif self.metaTranslator[pos]['typ'] == 'trunc1':
                        self.metaquery[self.metaTranslator[pos]['dst']] = products[key][pos][0:1]
                    else:
                        self.metaquery[self.metaTranslator[pos]['dst']] = products[key][pos]
                if pos in self.metaTranslator  and self.metaTranslator[pos]['tab']  in ['tiles','metatiles']:
                    self.tilequery[self.metaTranslator[pos]['dst']] = products[key][pos]
            #Get acquistion date from the filename
            acqdate = products[key]['beginposition'].date()
            acqtime = products[key]['beginposition'].time()
            self.acqdateStr = mj_dt.DateToStrDate(acqdate)
            #Get the day of year (doy) of acquistion 
            acqdoyStr = mj_dt.DateToYYYYDOY(acqdate)
            doy = int(acqdoyStr[4:7])
            
            self.tilequery['acqdate'] = acqdate
            self.tilequery['acqtime'] = acqtime
            self.tilequery['doy'] = doy
            
            self.tilequery['folder'] = self.tilequery['product']
            self.tilequery['source'] = products[key]['instrumentshortname']
            if self.tilequery['source'] == 'SAR-C SAR':
                self.tilequery['source'] = 'SAR-C'
                
            self._IdentifyMGRS(products,key)       

        #Register the search as done
        self.session._InsertVectorSearch(searchD)
           
    def _IdentifyMGRS(self, products, key):
        '''Identifies the MGRS
        '''
        #Get the 5th item from the identifier
        mgrsfn = products[key]['identifier'].split('_')[5]
        
        if self.metaquery['platformname'] == 'Sentinel-2' and 'tileid' in products[key]:
            print ('    Tile Id given', products[key]['tileid'])
            self.SetTile(products[key]['tileid'],products,key)

        elif self.metaquery['platformname'] == 'Sentinel-2' and len(mgrsfn) == 6 and mgrsfn[0] == 'T' and mgrsfn[1:2].isdigit():
            print ('    Tile Id from identifier', mgrsfn[1:6])
            self.SetTile(mgrsfn[1:6],products,key) 
            
        elif self.metaquery['platformname'] == 'Sentinel-1': 
            self._SetSentinelGranule(products,key)
        elif self.metaquery['platformname'] == 'Sentinel-2': 
            self._SetSentinelGranule(products,key)
        else:
            exit('error in sentine.sentinel _IdentifyMGRS')
                                      
    def _CreateGranulePolygon(self,products,key,mgrs):
        '''Create a shape vector (polygon) file 
        '''
        #The productgeom is 4236 (lon-lat)
        prodproj = mj_gis.MjProj()
        prodproj.SetFromEPSG(4326)
        productGeom = mj_gis.Geometry()
        productGeom.LoadWKTGeom(products[key]['footprint'])
        fieldDefD = {'type':'string','transfer':'constant','source':'globe','width':8}
        fieldDefL = [mj_gis.FieldDef('name',fieldDefD)]
        tarFP = os.path.join('/volumes',self.process.dstpath.volume,'sentinel','products','granules',products[key]['producttype'],self.acqdateStr)
        tarFN = '%(bn)s.shp' %{'bn':products[key]['identifier']}
        if not os.path.exists(tarFP):
            os.makedirs(tarFP)
        tarFPN = os.path.join(tarFP,tarFN)
        mj_gis.CreateESRIPolygonGeom(tarFPN, fieldDefL, productGeom, prodproj.proj_cs, mgrs)
        return productGeom       
      
    def SetTile(self, mgrs, products, key):
        utm = mgrs[0:2]
        mgrsid= mgrs[2:5]
        orbitid = products[key]['relativeorbitnumber']
        #Add the derived data to the tilequery dictionary
        self.tilequery['mgrs'] = mgrs
        self.tilequery['orbitid'] = orbitid
        self.tilequery['utm'] = utm
        self.tilequery['mgrsid'] = mgrsid       
        #Insert the meta and tile data in the postgres db
        self.session._InstertTileMeta(self.metaquery)
        self.session._InstertTile(self.tilequery)
        
    def _SetSentinelGranule(self, products, key):
        '''
        '''
        #Not an mgrs tile, this is a granule
        #Transsfer tilequery dict to granulequaery dict

        if self.metaquery['platformname'] == 'Sentinel-1':
            overlapTH = 0.01
        else:
            overlapTH = 0.5

        paramL = ['orbitid','acqdate','acqtime', 'doy', 'source', 'product', 'folder', 'filetype', 'filename']
        granulequery = {}
        granulequery['granuleid'] = products[key]['identifier']
        for item in self.tilequery:
            if item in paramL:
                granulequery[item] = self.tilequery[item]
        #Register the granule
        self.session._InstertGranule(granulequery)
        
        #Create a polygon and save as shape file
        granuleGeom = self._CreateGranulePolygon(products, key, 'granule')

        overlapD, mgrsBoundsD = self._OverlapGranuleTile(granuleGeom,overlapTH)

        #change tileid to granuleid as this is not a tile
        self.metaquery['granuleid'] = self.metaquery['tileid']
        self.metaquery.pop('tileid')
        #insert the granule
        self.session._InsertGranuleMeta(self.metaquery) 
        #Insert the MGRS tiles covered by this granule
        self.session._InsertGranuleTiles(granulequery['granuleid'], overlapD) 
      
    def _CheckPreviousSearches(self,startdate,enddate,platform,product,maxCloud):
        '''
        '''

        searchD = {}
        #minCloud = 0
        searchD['vectorfile'] = os.path.split(self.srcFPN)[1] 
        searchD['searchid'] = self.searchid
        searchD['platformname'] = platform
        searchD['product'] = product
        paramL = ['vectorfile', 'searchid', 'startdate','enddate','platformname','product','cloudcover']
        rec = self.session._SelectVectorSearch(searchD, paramL)

        if rec == None or self.process.overwrite:
            searchD['startdate'] = startdate
            searchD['enddate'] = enddate
            searchD['cloudcover'] = maxCloud
        else:            
            oldstartdate = rec[2]
            oldenddate = rec[3]
            oldcloudcover = rec[6]

            if oldstartdate <= startdate:
                startdate = oldenddate
            searchD['startdate'] = startdate

            if oldenddate >= enddate:
                enddate = oldstartdate

            searchD['enddate'] = enddate
            if oldcloudcover > maxCloud:
                maxCloud = oldcloudcover
            searchD['cloudcover'] = oldcloudcover
 
            if startdate > enddate and oldcloudcover >= maxCloud:
                #nothing to search for compared to earlier searhes
                return False
            elif startdate > enddate:
                print (searchD)
                BALLE
            else:
                return searchD
   
        #return startdate,enddate,platform,product,0,maxCloud,searchD
        return searchD

    def _ConstructTileLayer(self,tile):
        #TGTODO CHANGE TO DICT ISNTEAD of LIST
        #uuid, tileid, source, product, folder, utm, mgrsid, orbitid, acqdate = tile

        uuid, tileid, source, product, folder, acqdate, orbitid, utm, mgrsid, mgrs = tile
        #construct the composition
        compD = {'source':source, 'product':product, 'folder':folder, 'system':'sentinel', 'division':'tiles'}
        #Invoke the composition
        comp = SentinelComposition(compD)
        #Set the locus
        #locusD = {'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid}
        #Set the datum
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
        #Set the filename
        FN = '%(t)s.zip' %{'t':tileid}
        #Set the locus         
        loc = '%(utm)s%(mgrsid)s%(orbitid)s' %{'utm':utm,'mgrsid':mgrsid,'orbitid':orbitid}
        #Set the locuspath
        locusPath = os.path.join(utm,mgrsid)
        #Construct the locus dictionary
        locusD = {'locus':loc, 'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid, 'path':locusPath}
        #Invoke and return a SentinelTile             
        return SentinelTile(comp, locusD, datumD, self.process.dstpath, FN)
    
    def _ConstructGranuleLayer(self,granule):
        #uuid, tileid, source, product, folder, utm, mgrsid, orbitid, acqdate = tile
        print (granule)
        uuid, granuleid, source, product, folder, acqdate, orbitid = granule
        #construct the composition
        compD = {'source':source, 'product':product, 'folder':folder, 'system':'sentinel', 'division':'granules'}
        #Invoke the composition
        comp = SentinelComposition(compD)
        #Set the locus
        #locusD = {'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid}        #Set the datum
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
        #Set the filename
        FN = '%(t)s.zip' %{'t':granuleid}
        #Set the locus         
        #loc = '%(utm)s%(mgrsid)s%(orbitid)s' %{'utm':utm,'mgrsid':mgrsid,'orbitid':orbitid}
        #Set the locuspath
        locusPath = 'download'
        #Construct the locus dictionary
        locusD = {'locus':granuleid, 'region':granuleid, 'orbitid':orbitid, 'path':locusPath}
        #Invoke and return a SentinelTile             
        return SentinelTile(comp, locusD, datumD, self.process.dstpath, FN)

    def _GetGranuleTiles(self,granule):
        tarFP = os.path.join('/volumes',self.process.dstpath.volume,'sentinel','products','granules',prodkey['producttype'],acqdateStr)
        tarFN = '%(bn)s.shp' %{'bn':prodkey['identifier']}
        if not os.path.exists(tarFP):
            os.makedirs(tarFP)
        tarFPN = os.path.join(tarFP,tarFN)
        #print (tarFPN)

        mj_gis.CreateESRIPolygonGeom(tarFPN, fieldDefL, productGeom, prodproj.proj_cs, 'globe')
        print ('created product shape file',tarFPN)

    def _DownloadSentinelData(self):
        '''Download sentinel tiles and granules
        '''
        #Connect to the server
        api = SentinelAPI('thomas.gumbricht', 'iHg-98G-gug-34t', 'https://scihub.copernicus.eu/dhus')
        statusD = {}
        # TGTODO downloaded must be in xml, defaulted to N and not obligatory
        statusD['downloaded'] = self.process.params.downloaded
        #if self.process.params.platformname == 'Sentinel-1':  

        if self.process.proc.userProj.defregion == 'global':
            if self.process.params.orbitdirection.upper() != 'B':
                BALLE
            if self.process.params.tiles:
                tileL = self.session._SelectSentinelTiles(self.process.params,self.process.srcperiod,statusD)
                for tile in tileL:
                    print (tile)
                SNULLE
                self._GetTiles(tileL,api)
            else:
                granuleL = self.session._SelectSentinelGranules(self.process.params,self.process.srcperiod,statusD)
                self._GetGranules(granuleL,api)
                BALLE
        else:
            #print (self.process.proc.userProj.defregion)
            statusD['r.regionid'] = self.process.proc.userProj.defregion
            if self.process.params.orbitdirection.upper() != 'B':
                BALLE
            if self.process.params.tiles:
                
                tileL = self.session._SelectSentinelTiles(self.process.params,self.process.srcperiod,statusD)

                self._GetTiles(tileL,api)
            else:
                granuleL = self.session._SelectSentinelGranules(self.process.params,self.process.srcperiod,statusD)
                print ('granuleL',granuleL)
                SNULLE
                self._GetGranules(granuleL,api)
                
                BALLE
            BALLE

    def _DownloadSentinelTile(self,mgrs):
        '''
        '''
        #Connect to the server
        api = SentinelAPI('thomas.gumbricht', 'iHg-98G-gug-34t', 'https://scihub.copernicus.eu/dhus')
        statusD = {}
        # TGTODO downloaded must be in xml, defaulted to N and not obligatory
        statusD['downloaded'] = self.process.params.downloaded
        #if self.process.params.platformname == 'Sentinel-1':
        statusD['mgrs'] = mgrs
        print (statusD)
        tileL = self.session._SelectSentinelTiles(self.process.params,self.process.srcperiod,statusD)
        print (tileL)
        self._GetTiles(tileL,api)
            
    def _GetTiles(self,tileL,api):
        for tile in tileL:
            senTile = self._ConstructTileLayer(tile)
            uuid, tileid = tile[0:2]
            #Check if the tile exists
            senTile._Exists()
            if senTile._Exists():
                statusD = {'tileid': tileid,'column':'downloaded', 'status': 'Y'}
                self.session._UpdateTileStatus(statusD)
                statusD = {'tileid': tileid,'column':'organized', 'status': 'Y'}
                self.session._UpdateTileStatus(statusD)
                printstr = '    already downloaded %(d)s' %{'d':senTile.FN}
                print (printstr)
            else:
                printstr = '    downloading %(fpn)s' %{'fpn':senTile.FPN}
                print (printstr)

                api.download(uuid,senTile.FP)
                if os.path.exists(senTile.FPN):
                    statusD = {'tileid': tileid,'column':'downloaded', 'status': 'Y'}
                    self.session._UpdateTileStatus(statusD)
                    statusD = {'tileid': tileid,'column':'organized', 'status': 'Y'}
                    self.session._UpdateTileStatus(statusD)
                    
    def _GetGranules(self,granuleL,api):
        for granule in granuleL:
            senGranule = self._ConstructGranuleLayer(granule)
            uuid, granuleid = granule[0:2]
            #Check if the tile exists
            #senGranule._Exists()
            #print ('senGranule',senGranule._Exists() )
            #print ('senGranule', senGranule.FPN)
            #print ('FP', senGranule.FP)
            if senGranule._Exists():
                statusD = {'granuleid': granuleid,'column':'downloaded', 'status': 'Y'}
                self.session._UpdateGranuleStatus(statusD)
                printstr = '    already downloaded %(d)s' %{'d':senGranule.FN}
                print (printstr)
                #kolla om filen är duplicate
                for root, directories, filenames in os.walk('/Volumes/karttur2tb/sentinel/SAR-C/tiles/GRD'):
                    for filename in filenames:
                        if filename == senGranule.FN:
                            srcFPN = os.path.join(root,filename)
                            DELETECOPY
            else:
                printstr = '    downloading %(fpn)s' %{'fpn':senGranule.FPN}
                print ('senGranule.FN',senGranule.FN)
                downFlag = True
                #Check that the zip file did not end up under tiles:
                for root, directories, filenames in os.walk('/Volumes/karttur2tb/sentinel/SAR-C/tiles/GRD'):
                    for filename in filenames:
                        if filename == senGranule.FN:
                            srcFPN = os.path.join(root,filename)
                            downFlag = False
                            #move
                            printstr = 'Moving granule: %(src)s \n    to %(dst)s' %{'src': srcFPN, 'dst':senGranule.FPN}
                            print (printstr)
                            move(srcFPN,senGranule.FPN)
                            BALLE
    
                if downFlag:
                    return
                    api.download(uuid,senGranule.FP)
                    if os.path.exists(senGranule.FPN):
                        statusD = {'granuleid': granuleid,'column':'downloaded', 'status': 'Y'}
                        self.session._UpdateGranuleStatus(statusD)
                    else:
                        exit('Something fishy in _GetGranules')
                    
    def _ExplodeSentinel(self): 
        '''Explode downloaded sentinel data
        '''
 
        if self.process.params.platformname.lower() == 'sentinel-1': 
            self._ExplodeSentinel1Granules()
        elif self.process.params.platformname.lower() == 'sentinel-2':
            #Sentinel 2 can be both as tiles and as granules
            self._ExplodeSentinel2Tiles() 
            self._ExplodeSentinel2Granules()      
                   
    def _ExplodeSentinel2Granules(self):
        '''
        '''
        statusD = {}
        statusD['downloaded'] = 'Y'
        statusD['exploded'] = self.process.params.exploded
        if self.process.proj.defregion == 'global':
            if self.process.params.orbitdirection.upper() != 'B':
                NOTYET
            granuleL = self.session._SelectSentinelGranules(self.process.params,self.process.srcperiod,statusD)
        else:
            NOTYEAT
        for granule in granuleL:
            #construct the sentinel granule
            senGranule = self._ConstructGranuleLayer(granule)
            uuid, granuleid, source, product, folder, acqdate, orbitid = granule

            if not os.path.exists(senGranule.FPN):
                exitstr = 'EXTING, non-existing granule in _ExplodeSentinel2Granules'
                exit(exitstr)
            #Get the metatranslator to determine which metadata to save
            #self.metaTranslator = self.session._GetMetaTranslator()
            #self.metaL = [self.metaTranslator[key]['dst'] for key in self.metaTranslator if self.metaTranslator[key]['tab'] in ['meta','metatiles']]
            
            #Get all the meta for this granule
            metaRec = self.session._GetGranuleMeta(granuleid)

            if metaRec:
                paramD = ['product', 'proclevel', 'orbitnr', 'orbitdir', 'cloudcover', 'sensopmode', 's2datatakeid', 'procbase', 'platformid', 'platformname', 'instrument']
                self.metaquery = dict(zip(paramD,metaRec))
            else:
                exit('No metadata for granule in _ExplodeSentinel2Granules')
                
            granuleRec = self.session._GetGranuleTile(granuleid)
            if granuleRec:                
                paramD = ['orbitid', 'acqdate', 'acqtime', 'sunazimuth', 'sunelevation', 'doy', 'source', 'product', 'folder', 'filetype', 'filename', 'downloaded', 'organized', 'exploded', 'deleted', 'declouded', 'maskstatus', 'metacheck', 'tgnote']
                self.tilequery = dict(zip(paramD,granuleRec))
                popL = []
                for key in self.tilequery:
                    if self.tilequery[key] == None:
                        popL.append(key)
                for key in popL:
                    self.tilequery.pop(key)
            else:
                exit('No metadata for granule in _ExplodeSentinel2Granules')   

            granuleidPartsL = granuleid.split('_')

            sensingTime0 = granuleidPartsL[7]
            if sensingTime0[0] == 'V':
                sensingTime0 = sensingTime0[1:len(sensingTime0)]
            else:
                print ('sensingTime0',sensingTime0)
                exit ('Error in retrieving granule sensing stime')
            prodTime = granuleidPartsL[5]
            YYYYMMDD = mj_dt.DateToStrDate(self.tilequery['acqdate'])
            if YYYYMMDD != sensingTime0[0:8]:
                print ('YYYYMMDD != sensingTime0[0:8]',YYYYMMDD, sensingTime0[0:8])
                print ('YYYYMMDD',YYYYMMDD)
                print ('sensingTime0[0:8]',sensingTime0[0:8])
                exit ('Error in retrieving granule sensing stime')

            if product == 'S2MSI1C': 
                granulePath = self._ExplodeGranule(senGranule.FPN,'sentinel-2')
                statusD = {'granuleid': granuleid,'column':'exploded', 'status': 'Y'}
                self.session._UpdateGranuleStatus(statusD)
                #Get all granules
                tileL = os.listdir(granulePath)
                for tile in tileL:
                    #Skip any hidden files
                    if tile[0] == '.':
                        continue
                    self._ReadTileXml(granulePath,tile)

                    self.tilequery = dict(zip(paramD,granuleRec))
                    self.tilequery['sunazimuth'] = self.sunazimuth
                    self.tilequery['sunelevation'] = self.sunelevation
                    self.metaquery['cloudcover'] = self.cloudcover
                    tileid = tile
                    tileParts = tile.split('_')
                    mgrsfn = tileParts[len(tileParts)-2]
                    if len(mgrsfn) == 6 and mgrsfn[0] == 'T' and mgrsfn[1:2].isdigit():
                        print ('    Tile Id from identifier', mgrsfn[1:6])
                        mgrs = mgrsfn[1:6]
                        #MOVE TILE, REGISTER TILE AS DOWNLOADED AND EXPLODED
                    else:
                        exit('No mgrs identified for tile in _ExplodeSentinel2Granules')
            
                    procbase = tileParts[len(tileParts)-1]

                    if procbase[1:6] != self.metaquery['procbase']:
                        print ('procbase',procbase, procbase[1:6])
                        print ('self.metaquery',self.metaquery['procbase'])
                        exit('error in procbase')
                    else:
                        procbase = procbase.replace('.','')
                        
                    if self.metaquery['orbitnr'] < 10:
                        orbitnrStr = '00%(n)d' %{'n':self.metaquery['orbitnr']}
                    elif self.metaquery['orbitnr'] < 100:
                        orbitnrStr = '0%(n)d' %{'n':self.metaquery['orbitnr']}
                    else:
                        orbitnrStr = '%(n)d' %{'n':self.metaquery['orbitnr']}
                    
                    tileid = '%(MMM)s_MSIL1C_%(sensTime)s_%(procbase)s_R%(R000)s_%(Txxxxx)s_%(prodTime)s' %{'MMM':tileParts[0],
                                'sensTime':sensingTime0, 'procbase': procbase,'R000':orbitnrStr,'Txxxxx':mgrsfn,'prodTime':prodTime}
                    
                    self.metaquery['tileid'] = tileid
                    self.tilequery['tileid'] = tileid
                    
                    utm, mgrsid = mgrs[0:2], mgrs[2:5]
                    self.tilequery['utm'] = utm
                    self.tilequery['mgrsid'] = mgrsid
                    self.tilequery['mgrs'] = mgrs
                    
                    self.session._InstertTileMeta(self.metaquery)
                    self.session._InstertTile(self.tilequery)
                    
                    #construct the target, this is now an ordinary compositionprocess
                    compD = {'source':source,'product':product,'folder':'mask','band':'cloudmask','prefix':'cloudmask','suffix':'esa'}
                    comp = Composition(compD, self.process.system.srcsystem, self.process.system.srcdivision)
                    #Set the datum
                    datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
                    #Set the locus         
                    loc = '%(utm)s%(mgrsid)s' %{'utm':utm,'mgrsid':mgrsid}
                    #Set the locuspath
                    locusPath = os.path.join(utm,mgrsid)
                    #Construct the locus dictionary
                    locusD = {'locus':loc, 'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid, 'path':locusPath}

                    senTilePath = os.path.join(granulePath,tile)
                    #senTile = self._ConstructTileLayer(tile)
                    if not os.path.exists(senTilePath):
                        print ('senTilePath',senTilePath)
                        ERRORCRASH
                    #Register as exploded
                    statusD = {'tileid': tileid,'column':'exploded', 'status': 'Y'}
                    self.session._UpdateTileStatus(statusD)

                    self._OrganizeSentinel2GranuleTiles(senTilePath, comp, locusD, datumD, utm, mgrsid, acqdate, orbitid, tileid)
            
                statusD = {'granuleid': granuleid,'column':'organized', 'status': 'Y'}
                self.session._UpdateGranuleStatus(statusD)
            else:
                exitstr = 'Unknown product type in _ExplodeSentinel2Granules'
                exit(exitstr)
               
    def _ReadTileXml(self,granulePath,tile):
        xmlFP = os.path.join(granulePath,tile)
        #Get the xmlFN
        xmlFPL = os.listdir(xmlFP)
        for fitem in (xmlFPL):
            if fitem.endswith('.xml') and os.path.isfile(os.path.join(xmlFP,fitem)):
                tree = ET.parse(os.path.join(xmlFP,fitem))
                root = tree.getroot()
                for child in root:
                    if child.tag.split('}')[1] == 'Quality_Indicators_Info':
                        for grandchild in child:
                            if grandchild.tag == 'Image_Content_QI':
                                for greatgrandchild in grandchild:
                                    if greatgrandchild.tag == 'CLOUDY_PIXEL_PERCENTAGE':
                                        self.cloudcover = greatgrandchild.text
                    elif child.tag.split('}')[1] == 'Geometric_Info':
                        for grandchild in child:
                            if grandchild.tag == 'Tile_Angles':
                                for greatgrandchild in grandchild:
                                    if greatgrandchild.tag == 'Mean_Sun_Angle':
                                        for greatgreatgrandchild in greatgrandchild:
                                            if greatgreatgrandchild.tag == 'ZENITH_ANGLE':
                                                self.sunelevation = greatgreatgrandchild.text
                                            if greatgreatgrandchild.tag == 'AZIMUTH_ANGLE':
                                                self.sunazimuth = greatgreatgrandchild.text
           
    def _OrganizeSentinel2GranuleTiles(self, senTilePath, comp, locusD, datumD, utm, mgrsid, acqdate, orbitid, tileid):
        #Invoke the cloudmask as a standard vector layer  
        filepath = lambda: None
        filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = 'shp'
        cloudV = VectorLayer(comp, locusD, datumD, filepath)

        if self.process.params.setcloudmask:
            if not cloudV._Exists():
                compStr = 'MSK_CLOUDS'
                fileext = '.gml'
                cloudFPN = self._GetGMLMask(cloudV.FPN, senTilePath, compStr,fileext)
            else:
                cloudFPN = cloudV.FPN
        #"Construct dst path"
        queryD = {}
        queryD['product'] = {'val':self.process.params.prodtype, 'op':'=' }
        queryD['retrieve'] = {'val':'Y', 'op':'=' }
        extractL = self.session._SelectSentinelTemplate( queryD )
        for extract in extractL:
            source, product, folder, band, prefix, suffix, fileext, celltype, dataunit, resol, scalefac, offsetadd, cellnull, measure, retrieve, searchstr, bandname = extract 
            masked = 'Y'
            if not self.process.params.setcloudmask:
                folder = '%(f)s-cloudy' %{'f':folder}
                suffix = '%(s)s-cloudy' %{'s':suffix}
                masked = 'N'          
            compD = {'source':source,'product':product,'folder':folder,'band':band,'prefix':prefix,'suffix':suffix, 
                     'celltype':celltype, 'cellnull':cellnull, 'dataunit': folder, 'scalefac':scalefac, 'offsetadd':offsetadd, 'cellnull':cellnull, 'measure':measure,'masked':masked}
            comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
            datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
            
            #Set the locus         
            loc = '%(utm)s%(mgrsid)s' %{'utm':utm,'mgrsid':mgrsid}
            
            #Set the locuspath
            locusPath = os.path.join(utm,mgrsid)
            
            #Construct the locus dictionary
            locusD = {'locus':loc, 'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid, 'path':locusPath}
            
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = fileext
            
            #Create a standard reaster layer
            bandR = RasterLayer(comp, locusD, datumD, filepath)

            if bandR._Exists():
                print ('ALREADY DONE')
                self.session._InsertLayer(bandR)
                continue
  
            srcFPN  = self._GetBandFPN(senTilePath, searchstr,'.jp2')

            cmd = ['/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate']
            cmd.extend([ '-tr', '%(tr)d' %{'tr':resol} ,' %(tr)d' %{'tr':resol} ])
            cmd.extend(['-ot', celltype, '-a_nodata', '%(cn)d' %{'cn':cellnull} ])
            cmd.extend(['%(src)s' %{'src':srcFPN}, '%(dst)s' %{'dst':bandR.FPN} ]) 

            ThisProc = subprocess.check_call(cmd)
            print ('subprocess result', ThisProc)
            
            
            #Fix the cloud mask
            if self.process.params.setcloudmask and cloudFPN:
                print ('cloudFPN',cloudFPN)
                cloudLayer = os.path.splitext(os.path.split(cloudFPN)[1])[0]
                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_rasterize -burn %(cn)d -l %(l)s %(shp)s %(dst)s' %{'l':cloudLayer,'cn': cellnull, 'shp':cloudFPN, 'dst':bandR.FPN}
                os.system(cmd)
  
            #Fix all other masks
            maskBoolL = [self.process.params.setnodatamask, self.process.params.setdetfoomask, self.process.params.setdefectask, self.process.params.setsaturamask, self.process.params.settecquamask]
            maskIdL =  ['MSK_NODATA', 'MSK_DETFOO', 'MSK_DEFECT','MSK_SATURA', 'MSK_TECQUA']
            maskInvertL = [0,1,0,0,0]
            maskStatus = False
            for j, do in enumerate(maskBoolL):
                if do:  
                    maskStatus = True 
                    compStr = maskIdL[j]
                    endswithStr = '%(b)s_MSIL1C.gml' %{'b':bandname}

                    FN = '%(maskid)s_%(b)s.shp' %{'maskid':maskIdL[j],'b':bandname}
                    layerN = os.path.splitext(FN)[0]

                    dstFPN = os.path.join(os.path.split(cloudV.FPN)[0],FN)

                    if maskInvertL[j]:
                        
                        maskFPN = self._GetGMLMask(dstFPN, senTilePath, compStr, endswithStr, True)

                    else:
                        maskFPN = self._GetGMLMask(dstFPN, senTilePath, compStr, endswithStr)

                    if maskFPN and os.path.exists(maskFPN):
                        if maskIdL[j] not in ['MSK_NODATA','MSK_DETFOO','MSK_DEFECT']:
                            print ('nd',maskIdL[j])
                            exit("Error in _ExplodeS2MSI1C")
                        if maskInvertL[j]:
                            cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_rasterize -i -burn %(cn)d -l %(l)s %(shp)s %(dst)s' %{'l':layerN, 'cn': cellnull,'shp':maskFPN, 'dst':bandR.FPN}
                        else:
                            cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_rasterize -burn %(cn)d -l %(l)s %(shp)s %(dst)s' %{'l':layerN,'cn': cellnull, 'shp':maskFPN, 'dst':bandR.FPN}
                        os.system(cmd)
                        
            reclass = {}
            reclass[0] = {'op':'=', 'val': cellnull}
            replace = True
            RasterProcess.Reclass(bandR, bandR, reclass, replace) 
            self.session._InsertLayer(bandR)
            
        #Register as organized
        statusD = {'tileid': tileid,'column':'organized', 'status': 'Y'}
        self.session._UpdateTileStatus(statusD)
        #Register as declouded
        if self.process.params.setcloudmask:
            statusD = {'tileid': tileid,'column':'maskstatus', 'status': 'Y'}
            self.session._UpdateTileStatus(statusD)
            '''
            if maskStatus:
                statusD = {'tileid': tileid,'column':'maskstatus', 'status': 'Y'}
                self.session._UpdateTileStatus(statusD)
            '''
                     
    def _ExplodeGranule(self,zipFPN, platform):
        '''Explode the complete granule
        '''
        safePath =  zipper.ExplodeCompleteGranuleZip(zipFPN)
        if platform == 'sentinel-1':
            #granulePath = safePath
            granulePath = os.path.join(safePath,'measurement')
            '''
            tempPathL = os.listdir(tempPath)
            print (tempPathL)
            if len(tempPathL) == 1 and os.path.splitext(tempPathL[0])[1] == '.SAFE':
                if platform == 'sentinel-1':
                    granulePath = os.path.join(tempPath,tempPathL[0])
                elif platform == 'sentinel-2':
                    granulePath = os.path.join(tempPath,tempPathL[0],'GRANULE')
                else:
                    exit('undefined platform in _ExplodeGranule')
            else:
                exit('Can not find exploded path in _ExplodeGranule')
            '''
        else:
            granulePath = os.path.join(safePath,'GRANULE')
            '''
            subpath = os.path.split(zipFPN)[1]
            subpath = os.path.splitext(subpath)[0]
            subpath = '%(p)s.SAFE' %{'p':subpath}
            granulePath = os.path.join(tempPath,subpath)
            print ('granulepath',granulePath)
            BALLE
            '''
        if os.path.isdir(granulePath):
            return granulePath
        else:
            exit('granulePath does not exist in _ExplodeGranule')
                          
    def _ExplodeSentinel2Tiles(self):
        '''Explode sentinel tiles
        '''
        statusD = {}
        statusD['downloaded'] = 'Y'
        statusD['exploded'] = self.process.params.exploded
        if self.process.proj.defregion == 'global':
            if self.process.params.orbitdirection.upper() != 'B':
                BALLE
            tileL = self.session._SelectSentinelTiles(self.process.params,self.process.srcperiod,statusD)
        else:
            BALLE
        for tile in tileL:
            uuid, tileid, source, product, folder, acqdate, orbitid, utm, mgrsid, mgrs  = tile

            if uuid == None:
                #uuid == None is true for tiles that are retrieved from granules, and they are handled differencely
                pass
            else:
                #construct the sentinel tile
                senTile = self._ConstructTileLayer(tile)
                if not os.path.exists(senTile.FPN):
                    print (senTile.FPN)
                    MISSINGTILE
               
                #construct the target, this is now an ordinary compositionprocess
                compD = {'source':source,'product':product,'folder':'mask','band':'cloudmask','prefix':'cloudmask','suffix':'esa'}
                comp = Composition(compD, self.process.system.srcsystem, self.process.system.srcdivision)
                #Set the datum
                datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
                #Set the locus         
                loc = '%(utm)s%(mgrsid)s' %{'utm':utm,'mgrsid':mgrsid}
                #Set the locuspath
                locusPath = os.path.join(utm,mgrsid)
                #Construct the locus dictionary
                locusD = {'locus':loc, 'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid, 'path':locusPath}
                print ('product',product)
    
                #Processing depends on data product
                if product == 'S2MSI1C':  
                    self._ExplodeS2MSI1C(senTile, comp, locusD, datumD, utm, mgrsid, acqdate, orbitid, tileid)
                else:
                    exitstr = 'Unknown product type in _ExplodeSentinel2Tiles'
                    exit(exitstr)
  
    def _ExplodeS2MSI1C(self, senTile, comp, locusD, datumD, utm, mgrsid, acqdate, orbitid, tileid):
        #Invoke the cloudmask as a standard vector layer  
        filepath = lambda: None
        filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = 'shp'
        cloudV = VectorLayer(comp, locusD, datumD, filepath)
        if self.process.params.setcloudmask:
            if not cloudV._Exists():
                compStr = 'MSK_CLOUDS_B00.gml'
                cloudFPN = self._ExtractGMLMask(cloudV.FPN, senTile.FPN, compStr)
            else:
                cloudFPN = cloudV.FPN
        #"Construct dst path"
        queryD = {}
        queryD['product'] = {'val':self.process.params.prodtype, 'op':'=' }
        queryD['retrieve'] = {'val':'Y', 'op':'=' }
        extractL = self.session._SelectSentinelTemplate( queryD )
        for extract in extractL:
            source, product, folder, band, prefix, suffix, fileext, celltype, dataunit, resol, scalefac, offsetadd, cellnull, measure, retrieve, searchstr, bandname = extract 
            masked = 'Y'
            if not self.process.params.setcloudmask:
                folder = '%(f)s-cloudy' %{'f':folder}
                suffix = '%(s)s-cloudy' %{'s':suffix}
                masked = 'N'          
            compD = {'source':source,'product':product,'folder':folder,'band':band,'prefix':prefix,'suffix':suffix, 
                     'celltype':celltype, 'cellnull':cellnull, 'dataunit': folder, 'scalefac':scalefac, 'offsetadd':offsetadd, 'cellnull':cellnull, 'measure':measure,'masked':masked}
            comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
            datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
            
            #Set the locus         
            loc = '%(utm)s%(mgrsid)s' %{'utm':utm,'mgrsid':mgrsid}
            
            #Set the locuspath
            locusPath = os.path.join(utm,mgrsid)
            
            #Construct the locus dictionary
            locusD = {'locus':loc, 'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid, 'path':locusPath}
            
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = fileext
            
            #Create a standard reaster layer
            bandR = RasterLayer(comp, locusD, datumD, filepath)

            if bandR._Exists():
                print ('ALREADY DONE')
                self.session._InsertLayer(bandR)
                continue
            
            srcFPN  = zipper.UnZip(senTile.FPN, searchstr)
            cmd = ['/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate']
            cmd.extend([ '-tr', '%(tr)d' %{'tr':resol} ,' %(tr)d' %{'tr':resol} ])
            cmd.extend(['-ot', celltype, '-a_nodata', '%(cn)d' %{'cn':cellnull} ])
            cmd.extend(['%(src)s' %{'src':srcFPN}, '%(dst)s' %{'dst':bandR.FPN} ]) 

            ThisProc = subprocess.check_call(cmd)
            print ('subprocess result', ThisProc)
            
            #Fix the cloud mask
            if self.process.params.setcloudmask and cloudFPN:
                cloudLayer = os.path.splitext(os.path.split(cloudFPN)[1])[0]
                cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_rasterize -burn %(cn)d -l %(l)s %(shp)s %(dst)s' %{'l':cloudLayer,'cn': cellnull, 'shp':cloudFPN, 'dst':bandR.FPN}
                os.system(cmd)
                
            #Fix all other masks
            maskBoolL = [self.process.params.setnodatamask, self.process.params.setdetfoomask, self.process.params.setdefectask, self.process.params.setsaturamask, self.process.params.settecquamask]
            maskIdL =  ['MSK_NODATA', 'MSK_DETFOO', 'MSK_DEFECT','MSK_SATURA', 'MSK_TECQUA']
            maskInvertL = [0,1,0,0,0]
            for j, do in enumerate(maskBoolL):
                if do:
                    compStr = '%(maskid)s_%(b)s.gml' %{'maskid':maskIdL[j],'b':bandname}
                    layerN = os.path.splitext(compStr)[0]

                    FN = compStr.replace('.gml','.shp')
                    
                    dstFPN = os.path.join(os.path.split(cloudV.FPN)[0],FN)

                    if maskInvertL[j]:
                        maskFPN = self._ExtractGMLMask(dstFPN, senTile.FPN, compStr, True)
                    else:
                        maskFPN = self._ExtractGMLMask(dstFPN, senTile.FPN, compStr)

                    if maskFPN and os.path.exists(maskFPN):
                        if maskIdL[j] not in ['MSK_NODATA','MSK_DETFOO','MSK_DEFECT']:
                            print ('nd',maskIdL[j])
                            exit("Error in _ExplodeS2MSI1C")
                        if maskInvertL[j]:
                            cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_rasterize -i -burn %(cn)d -l %(l)s %(shp)s %(dst)s' %{'l':layerN, 'cn': cellnull,'shp':maskFPN, 'dst':bandR.FPN}
                        else:
                            cmd = '/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_rasterize -burn %(cn)d -l %(l)s %(shp)s %(dst)s' %{'l':layerN,'cn': cellnull, 'shp':maskFPN, 'dst':bandR.FPN}
                        os.system(cmd)
                        
            reclass = {}
            reclass[0] = {'op':'=', 'val': cellnull}
            RasterProcess.Reclass(bandR, bandR, reclass) 
            self.session._InsertLayer(bandR)
            
        #Register as exploded
        statusD = {'tileid': tileid,'column':'exploded', 'status': 'Y'}
        self.session._UpdateTileStatus(statusD)
         
    def _ExplodeSentinel1Granules(self):
        overlapTH = self.process.params.granuleoverlap
        statusD = {}
        statusD['downloaded'] = 'Y'
        #statusD['exploded'] = self.process.params.exploded
        if self.process.proj.defregion == 'global':
            if self.process.params.orbitdirection.upper() != 'B':
                NOTYET
            granuleL = self.session._SelectSentinelGranules(self.process.params,self.process.srcperiod,statusD)
        else:
            NOTYEAT
        for granule in granuleL:
            #construct the sentinel granule
            senGranule = self._ConstructGranuleLayer(granule)
            uuid, granuleid, source, product, folder, acqdate, orbitid = granule

            if not os.path.exists(senGranule.FPN):
                exitstr = 'EXTING, non-existing granule in _ExplodeSentinel2Granules'
                exit(exitstr)
            #Get the metatranslator to determine which metadata to save
            #self.metaTranslator = self.session._GetMetaTranslator()
            #self.metaL = [self.metaTranslator[key]['dst'] for key in self.metaTranslator if self.metaTranslator[key]['tab'] in ['meta','metatiles']]
            
            #Get all the meta for this granule
            metaRec = self.session._GetGranuleMeta(granuleid)
            if metaRec:
                paramD = ['product', 'proclevel', 'orbitnr', 'orbitdir', 'cloudcover', 'sensopmode', 's2datatakeid', 'procbase', 'platformid', 'platformname', 'instrument']
                self.metaquery = dict(zip(paramD,metaRec))
            else:
                exit('No metadata for granule in _ExplodeSentinel2Granules')
                
            granuleRec = self.session._GetGranuleTile(granuleid)
            print ('granuleRec',granuleRec)
            if granuleRec:                
                paramD = ['orbitid', 'acqdate', 'acqtime', 'sunazimuth', 'sunelevation', 'doy', 'source', 'product', 'folder', 'filetype', 'filename', 'downloaded', 'organized', 'exploded', 'deleted', 'declouded', 'maskstatus', 'metacheck', 'tgnote']
                self.tilequery = dict(zip(paramD,granuleRec))
                popL = []
                for key in self.tilequery:
                    if self.tilequery[key] == None:
                        popL.append(key)
                for key in popL:
                    self.tilequery.pop(key)
            else:
                exit('No metadata for granule in _ExplodeSentinel2Granules')   

            granuleidPartsL = granuleid.split('_')
            missionid, mode, prodtyp, prodclass, startdatetime, enddatetime, absorbnr, datatakeid, prodid = granuleidPartsL
            #print (granuleidPartsL)
            #print (missionid, mode, prodtyp, prodclass, startdatetime, enddatetime, absorbnr, datatakeid, prodid)
            resolclass = prodtyp[3]
            prodtyp = prodtyp[0:3]
            proclevel = prodclass[0]
            polarisation = prodclass[2:4]
            prodclass = prodclass[1]
            #print (resolclass, prodtyp, proclevel, polarisation, prodclass)

            YYYYMMDD = mj_dt.DateToStrDate(self.tilequery['acqdate'])
            if YYYYMMDD != startdatetime[0:8]:
                print ('YYYYMMDD != sensingTime0[0:8]',YYYYMMDD, startdatetime[0:8])
                print ('YYYYMMDD',YYYYMMDD)
                print ('sensingTime0[0:8]',startdatetime[0:8])
                exit ('Error in retrieving granule sensing stime')
            if product != prodtyp:
                print ('product != prodtyp',product, prodtyp)
                exit ('Error in retrieving sentinel 1 granule product')

            if product == 'GRD': 
                granulePath = self._ExplodeGranule(senGranule.FPN, 'sentinel-1')
                #statusD = {'granuleid': granuleid,'column':'exploded', 'status': 'Y'}
                #self.session._UpdateGranuleStatus(statusD)
                #Get the tiles to construct
                tileL = self.session._SelectGranuleTiles(granuleid, overlapTH)

                #Loop over the tiles to extract for this granule
                for tile in tileL:

                    #Get the mgrs details
                    rec = self.session._SelectMGRS(tile)
                    if rec == None:
                        #exit('No mgrs found')
                        pass
                    else:
                        utmzone, mgrsid, proj4, minx, miny, maxx, maxy, refsize, refcols, reflins = rec
                        #Now I have to get the bands
                        #print ('projection',utmzone, mgrsid, proj4, minx, miny, maxx, maxy, refsize, refcols, reflins)
                        #print ('granulepath',granulePath)
                        self._ExtractGRDTile(granulePath, acqdate, orbitid, utmzone, mgrsid, proj4, minx, miny, maxx, maxy, refsize, refcols, reflins)
            
                statusD = {'granuleid': granuleid,'column':'downloaded', 'status': 'Y'}
                self.session._UpdateGranuleStatus(statusD)
                statusD = {'granuleid': granuleid,'column':'exploded', 'status': 'Y'}
                self.session._UpdateGranuleStatus(statusD)
                statusD = {'granuleid': granuleid,'column':'organized', 'status': 'Y'}
                self.session._UpdateGranuleStatus(statusD)
            
    def _ExtractGRDTile(self, granuleFP, acqdate, orbitid, utmzone, mgrsid, proj4, minx, miny, maxx, maxy, refsize, refcols, reflins):
        #"Construct dst path"
        queryD = {}
        queryD['product'] = {'val':self.process.params.prodtype, 'op':'=' }
        queryD['retrieve'] = {'val':'Y', 'op':'=' }
        extractL = self.session._SelectSentinelTemplate( queryD )
        print ('extractL', extractL)
        if utmzone < 10:
            utm = '0%(utm)d' %{'utm':utmzone}
        else:
            utm = '%(utm)d' %{'utm':utmzone}
        for extract in extractL:
            #I have to get the exploded source file 
            print ('extract',extract)
            #get a list of files
            FN = False
            fileL = os.listdir(granuleFP)
            for file in fileL:
                if extract[15] in file:
                    FN = file
            if not FN:
                exit('No band file found in sentinel.sentinel _ExtractGRDTile')        
            granuleFPN = os.path.join(granuleFP, FN)

            source, product, folder, band, prefix, suffix, fileext, celltype, dataunit, resol, scalefac, offsetadd, cellnull, measure, retrieve, searchstr, bandname = extract 
            compD = {'source':source,'product':product,'folder':folder,'band':band,'prefix':prefix,'suffix':suffix, 
                     'celltype':celltype, 'cellnull':cellnull, 'dataunit': folder, 'scalefac':scalefac, 'offsetadd':offsetadd, 'cellnull':cellnull, 'measure':measure,'masked':'Y'}
            comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
            datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
            #Set the locus         
            #loc = '%(utm)s%(mgrsid)s%(orbitid)s' %{'utm':utm,'mgrsid':mgrsid,'orbitid':orbitid}

            loc = '%(utm)s%(mgrsid)s' %{'utm':utm,'mgrsid':mgrsid}
            #Set the locuspath
            
            locusPath = os.path.join(utm,mgrsid)
            #Construct the locus dictionary
            locusD = {'locus':loc, 'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid, 'path':locusPath}
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = fileext
            bandR = RasterLayer(comp, locusD, datumD, filepath)
            print ('    bandR.FPN',bandR.FPN)

            if bandR._Exists():
                self.session._InsertLayer(bandR)
                continue
            
            '''
            cmd = ['/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdal_translate']
            cmd.extend([ '-tr', '%(tr)d' %{'tr':resol} ,' %(tr)d' %{'tr':resol} ])
            cmd.extend(['-ot', celltype, '-a_nodata', '%(cn)d' %{'cn':cellnull} ])
            cmd.extend(['%(src)s' %{'src':granuleFPN}, '%(dst)s' %{'dst':bandR.FPN} ]) 
            print (cmd)
            '''
            cmd = ['/Library/Frameworks/GDAL.framework/Versions/2.1/Programs/gdalwarp']
            cmd.extend([ '-t_srs', '%(tsrs)s' %{'tsrs':proj4}])
            cmd.extend([ '-te', '%(xmin)d' %{'xmin':minx, 'ymin':miny, 'xmax':maxx, 'ymax':maxy} , '%(ymin)d' %{'xmin':minx, 'ymin':miny, 'xmax':maxx, 'ymax':maxy} ,'%(xmax)d' %{'xmin':minx, 'ymin':miny, 'xmax':maxx, 'ymax':maxy} , '%(ymax)d' %{'xmin':minx, 'ymin':miny, 'xmax':maxx, 'ymax':maxy} ])
            cmd.extend([ '-tr', '%(tr)d' %{'tr':resol} ,' %(tr)d' %{'tr':refsize} ])
            cmd.extend(['-ot', celltype,  '-srcnodata', '0', '-dstnodata', '%(cn)d' %{'cn':cellnull} ])
            cmd.extend(['%(src)s' %{'src':granuleFPN}, '%(dst)s' %{'dst':bandR.FPN} ]) 

            ThisProc = subprocess.check_call(cmd)
            self.session._InsertLayer(bandR)
     
    def _ExtractGMLMask(self, dstFPN, zipFPN, compStr, union = False):
        '''
        '''
        #def ExtractGMLMask(dstFP, zipFPN, compStr, pathId, rowId, datumStr, relorbit, union = False, folder = 'mask'):
        srcFPN  = zipper.UnZip(zipFPN, compStr)
        #maskFN = os.path.split(srcFPN)[1].replace('.gml','.shp')
        #maskFP = os.path.split(dstFP)[0]
        #print 'maskFP',mdstFPN
    
        #maskFP = os.path.join(maskFP,folder,pathId,rowId,datumStr,relorbit)
        #print 'maskFP',maskFP
        #maskFPN = os.path.join(maskFP, maskFN)
        #if not os.path.exists(maskFP):
        #    os.makedirs(maskFP)

        maskFPN = self._get_mask_shp(srcFPN, dstFPN, union)
        #print ('maskFPN',maskFPN)
        return maskFPN
    
    def _GetGMLMask(self, dstFPN, srcFP, compStr, fileext, union = False):
        '''
        '''
        #Get all the files under the topfolder
        for root, directories, filenames in os.walk(srcFP):
            for filename in filenames:
                print ('filename',filename, compStr, fileext) 
                if compStr in filename and filename.endswith(fileext):
                    srcFPN = os.path.join(root,filename)
                    print ( srcFPN )
                    maskFPN = self._get_mask_shp(srcFPN, dstFPN, union)
                    print ( dstFPN )
                    return maskFPN
        BALLE
        return False
        BALLE
        #maskFN = os.path.split(srcFPN)[1].replace('.gml','.shp')
        #maskFP = os.path.split(dstFP)[0]
        #print 'maskFP',mdstFPN
    
        #maskFP = os.path.join(maskFP,folder,pathId,rowId,datumStr,relorbit)
        #print 'maskFP',maskFP
        #maskFPN = os.path.join(maskFP, maskFN)
        #if not os.path.exists(maskFP):
        #    os.makedirs(maskFP)

        maskFPN = self._get_mask_shp(srcFPN, dstFPN, union)
        #print ('maskFPN',maskFPN)
        return maskFPN
    
    def _GetBandFPN(self, srcFP, compStr, fileext):
        '''
        '''
        #Get all the files under the topfolder
        for root, directories, filenames in os.walk(srcFP):
            for filename in filenames:
                print ('filename',filename) 
                if compStr in filename and os.path.splitext(filename)[1] == fileext:
                    srcFPN = os.path.join(root,filename)
                    print ( srcFPN )
                    return srcFPN 
        return False
        BALLE
        #maskFN = os.path.split(srcFPN)[1].replace('.gml','.shp')
        #maskFP = os.path.split(dstFP)[0]
        #print 'maskFP',mdstFPN
    
        #maskFP = os.path.join(maskFP,folder,pathId,rowId,datumStr,relorbit)
        #print 'maskFP',maskFP
        #maskFPN = os.path.join(maskFP, maskFN)
        #if not os.path.exists(maskFP):
        #    os.makedirs(maskFP)

        maskFPN = self._get_mask_shp(srcFPN, dstFPN, union)
        #print ('maskFPN',maskFPN)
        return maskFPN
    
    def _get_mask_shp(self, path_in_gml, path_out, union = False):
        """Transform the gml cloud mask to a shape file."""
               
        name = os.path.split(path_in_gml)[-1]
    
        if name.endswith('gml'):
            #path_out_shp = os.path.join(path_out,name.replace('.gml','.shp'))
                    
            src_srs = self._get_gml_src_proj(path_in_gml)
    
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
            
    def _get_gml_src_proj(self, path_gml):
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
    
    def _FindGranuleTiles(self):
        statusD = {}
        overlapTH = self.process.params.granuleoverlap

        statusD['downloaded'] = 'Y'
        statusD['exploded'] = 'Y'
        if self.process.proj.defregion == 'global':
            granuleL = self.session._SelectSentinelGranules(self.process.params,self.process.srcperiod,statusD)
        else:
            NOTYEAT
        for granule in granuleL:
            uuid, granuleid, source, product, folder, acqdate, orbitid = granule
            print ('granule',granule)
            self.acqdateStr = mj_dt.DateToStrDate(acqdate)
            shpFP = os.path.join('/volumes',self.process.dstpath.volume,'sentinel','products','granules',self.process.params.prodtype,self.acqdateStr)
            shpFN = '%(bn)s.shp' %{'bn':granuleid}

            shpFPN = os.path.join(shpFP,shpFN)
            if not os.path.isfile(shpFPN):   
                exitstr =  'In sentinel.sentinel missing granule shape file: %(s)s' %{'e': shpFPN}
                exit(exitstr)
            print (shpFPN)
            srcDS,srcLayer = mj_gis.ESRIOpenGetLayer(shpFPN)

            #Loop over the features in the layer
            for feature in srcLayer.layer:
                #Craete a geom instance
                granuleGeom = mj_gis.Geometry()
                #add the feature and extract the geom
                granuleGeom.GeomFromFeature(feature)
                if srcLayer.geomtype.lower() != 'polygon':
                    exit('Expected a polygon in _FindGranuleTiles')
            overlapD, mgrsBoundsD =  self._OverlapGranuleTile(granuleGeom,overlapTH)
            STPP #DO SOMETHING
        
    def _OverlapGranuleTile(self,granuleGeom,overlapTH):                 
        #Get the bounds of the granule geom
        west, south, east, north = granuleGeom.shapelyGeom.bounds
        #Find out which tiles are covered by this granule
        '''
        sentinelflag = 'Y' # 'Y' indicates that only tiles that are included in the sentinel system are retrieved
        recs2 = self.session._SearchMGRSFromWSEN(west, south, east, north, sentinelflag)
        print (len(recs2))
        recs2 = self.session._SearchMGRSFromWSEN(west, south, east, north, False)
        print (len(recs2))
        '''
        recs = self.session._SearchTilesFromWSEN(west, south, east, north)
        
        if len(recs) == 0:
            exit('No granule tile match in sentine.sentinel') 
        
        #Loop over the recs and find out which of them fits inside the product polygon
        overlapD = {}
        mgrsBoundsD = {}
        paramL = ['minx', 'miny', 'maxx', 'maxy']
        for tile in recs:  
            mgrs,tilewest,tilesouth,tileeast,tilenorth,ullon,ullat,urlon,urlat,lrlon,lrlat,lllon,lllat, minx, miny, maxx, maxy = tile  
           
            llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
            tilegeom = mj_gis.Geometry()
            tilegeom.PointsToPolygonGeom(llptL)
            #Get the overlap
            overlapGeom = tilegeom.ShapelyIntersection(granuleGeom)
            productoverlap = overlapGeom.area/tilegeom.shapelyGeom.area
            valueL = [minx, miny, maxx, maxy]
            if productoverlap >= overlapTH:
                overlapD[mgrs] = productoverlap
                mgrsBoundsD[mgrs] = dict(zip(paramL,valueL))
        return overlapD, mgrsBoundsD

    def _ReOrganiseSentinel(self):
        searchpath = self.process.params.searchpath
        searchtype = '.zip'
        for root, directories, filenames in os.walk(searchpath):
            for filename in filenames:
                if filename.endswith(searchtype) and os.path.isfile(os.path.join(root,filename)):
                    srcFPN = os.path.join(root,filename)
                    if self.process.params.tiles:
                        BALLE
                        tileL = self.session._SelectSentinelTiles(self.process.params,self.process.srcperiod,statusD)
                        self._GetTiles(tileL,api)
                    else:                        
                        granuleid = os.path.splitext(filename)[0]
                        statusD = {'granuleid':granuleid}
                        granuleL = self.session._SelectSentinelGranules(self.process.params,False,statusD)
                        print ('granuleL',granuleL)
                        if len(granuleL) == 1:
                            senGranule = self._ConstructGranuleLayer(granuleL[0])
                            if senGranule._Exists():
                                os.remove(srcFPN)

                            else:
                                #move
                                printstr = 'Moving granule: %(src)s \n    to %(dst)s' %{'src': srcFPN, 'dst':senGranule.FPN}
                                print (printstr)
                                move(srcFPN,senGranule.FPN)
                                statusD = {'granuleid': granuleid,'column':'downloaded', 'status': 'Y'}
                                self.session._UpdateGranuleStatus(statusD)

                        elif len(granuleL) == 0:
                            granuleidPartsL = granuleid.split('_')
                            print (granuleidPartsL)
                            BALLE
                        else:
                            BALLE


    def _GeoCheckSentinelTilesOld(self):
        '''Check the geometry of Sentinel tiles
        '''
        statusD = {}
        statusD['exploded'] = 'Y'
        if self.process.proj.defregion == 'global':
            if self.process.params.orbitdirection.upper() != 'B':
                NOTYET
            tileL = self.session._SelectSentinelTiles(self.process.params,self.process.srcperiod,statusD)
        else:
            NOTYET
        for tile in tileL:
            uuid, tileid, source, product, folder, acqdate, orbitid, utm, mgrsid, mgrs = tile
            #"Construct dst path"
            queryD = {}
            queryD['product'] = {'val':self.process.params.prodtype, 'op':'=' }
            queryD['retrieve'] = {'val':'Y', 'op':'=' }
            extractL = self.session._SelectSentinelTemplate( queryD )
            for extract in extractL:
                source, product, folder, band, prefix, suffix, fileext, celltype, dataunit, resol, scalefac, offsetadd, cellnull, measure, retrieve, searchstr, bandname = extract 
                compD = {'source':source,'product':product,'folder':folder,'band':band,'prefix':prefix,'suffix':suffix} 
                comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
                datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}
                #Set the locus         
                loc = '%(utm)s%(mgrsid)s' %{'utm':utm,'mgrsid':mgrsid}
                #Set the locuspath
                locusPath = os.path.join(utm,mgrsid)
                #Construct the locus dictionary
                locusD = {'locus':loc, 'utm':utm, 'mgrsid':mgrsid, 'orbitid':orbitid, 'path':locusPath}
                #Construct the filepath
                filepath = lambda: None
                filepath.volume = self.process.srcpath.volume; filepath.hdrfiletype = fileext
                bandR = RasterLayer(comp, locusD, datumD, filepath)
                if bandR._Exists():
                    queryD = {'mgrs':mgrs}
                    #print (bandR.FPN)
                    #Open and read layer metadata for lins,cols and projection
                    bandR.GetRastermetadata()
                    #print ('spatialRef',bandR.metadata.cols,bandR.metadata.lins, bandR.metadata.bounds)
                    #print ('metadata',bandR.metadata)
                    minx = bandR.metadata.bounds[0]
                    maxx = bandR.metadata.bounds[2]
                    miny = bandR.metadata.bounds[1]
                    maxy = bandR.metadata.bounds[3]
                    queryD['minx'] = minx
                    queryD['maxx'] = maxx
                    queryD['miny'] = miny
                    queryD['maxy'] = maxy
                    if bandR.metadata.cellsize == 20.0:
                        queryD['refcols'] = bandR.metadata.cols
                        queryD['reflins'] = bandR.metadata.lins    
                    else:
                        BALLE
                    #
                    ptL = ( (minx,maxy),(maxx,maxy),(maxx,miny),(minx,miny) )
                    cornergeom = mj_gis.Geometry()
                    cornergeom.PointsToMultiPointGeom(ptL)
                    #Reproject from UTM to lonlat to get the correct corners in lonlat
                    srcproj = mj_gis.MjProj()
                    srcproj.SetFromProj4(bandR.spatialRef.proj4)
                    dstproj = mj_gis.MjProj()
                    dstproj.SetFromEPSG(4326)
                    lonlatgeom = srcproj.ReprojectGeom(cornergeom, dstproj)
                    #ptL = ((query['ullon'],query['ullat']), (query['urlon'],query['urlat']), (query['lrlon'],query['lrlat']), ((query['lllon'],query['lllat'])) )
                    #Set ul and lr corners of bounding box

                    lonlat = [(p.x, p.y) for p in lonlatgeom.shapelyGeom]
                    queryD['ullat'] = lonlat[0][1]
                    queryD['ullon'] = lonlat[0][0]
                    queryD['urlat'] = lonlat[1][1]
                    queryD['urlon'] = lonlat[1][0]
                    queryD['lrlat'] = lonlat[2][1]
                    queryD['lrlon'] = lonlat[2][0]
                    queryD['lllat'] = lonlat[3][1]
                    queryD['lllon'] = lonlat[3][0]
                
                    queryD['cellsize'] = 20.0
                    queryD['projection'] = bandR.metadata.projection
                    queryD['verified'] = 'Y'  
                    print (queryD)
                    print ('spatialref',bandR.spatialRef.proj_cs)
                    print ('proj4',bandR.spatialRef.proj4)

                    
                    utmgeom = dstproj.ReprojectGeom(lonlatgeom, srcproj)
                    utm = [(p.x, p.y) for p in utmgeom.shapelyGeom]
                    print ('utm',utm)
                    
                    BALLE
                    self.session._InsertTileCoords(queryD)
                    '''
                    print ('lins',bandR.metadata.lins)
                    print ('cols',bandR.metadata.cols)
                    print ('cellsize',bandR.metadata.cellsize)
                    print ('bounds',bandR.metadata.bounds)
                    print ('gt',bandR.metadata.geotrans)
                    print ('projection',bandR.metadata.projection)

                    
                    self.lins = self.datasource.RasterYSize
                    self.cols = self.datasource.RasterXSize
                    self.cellnull = self.layer.GetNoDataValue()
                    self.celltype = gdal.GetDataTypeName(self.layer.DataType)
                    self.projection = self.datasource.GetProjection()
                    self.geotrans = self.datasource.GetGeoTransform()
            
                    #Get the extent of the image
                    self.ext=[]
                    xarr=[0,self.cols]
                    yarr=[0,self.lins]
                    for px in xarr:
                        for py in yarr:
                            x=self.gt[0]+(px*self.gt[1])+(py*self.gt[2])
                            y=self.gt[3]+(px*self.gt[4])+(py*self.gt[5])
                            self.ext.append([x,y])
                        yarr.reverse()
                    self.bounds = (self.ext[0][0], self.ext[2][1], self.ext[2][0],self.ext[0][1])
                    #Get the spatial resolution
                    cellsize = [(self.ext[2][0]-self.ext[0][0])/self.cols, (self.ext[0][1]-self.ext[2][1])/self.lins] 
                    if cellsize[0] != cellsize[1]:
                        pass
                    self.cellsize = cellsize[0]
                    '''
                    COOL
                
                
                
    def _LinkDefaultRegionsToSentinel(self):
        '''This is just a complicated manner to get the sentinel tile file in lonlat projection
        '''
        for locus in self.process.srcLayerD:
            if len(self.process.srcLayerD[locus]) == 0:
                exitstr = 'EXITING, no dates defined in Sentinel._ExtractSentinelTileCoords'
                exit(exitstr)
            for datum in self.process.srcLayerD[locus]:
                if len(self.process.srcLayerD[locus][datum]) == 0:
                    exitstr = 'EXITING, no compositions defined in Sentinel._ExtractSentinelTileCoords'
                    exit(exitstr)
                for comp in self.process.srcLayerD[locus][datum]:
                    self.srcLayer = self.process.srcLayerD[locus][datum][comp]
        self._GetSentinelTilesDict()
        self._GetSystemDefRegions()

        
    def _GetSentinelTilesDict(self):
        '''
        '''
        recs = self.session._SelectSentinelTileCoords({})
        self.senTileD ={}
        for rec in recs:
            epsg,mgrs,utmzone,mgrsid,minx,miny,maxx,maxy,ullat,ullon,lrlat,lrlon,urlat,urlon,lllat,lllon = rec
            llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
            sentilegeom = mj_gis.Geometry()
            sentilegeom.PointsToPolygonGeom(llptL)
            west, south, east, north = sentilegeom.shapelyGeom.bounds
            self.senTileD[mgrs] = {'mgrs':mgrs,'mgrsid':mgrsid,'utmzone':utmzone,'geom':sentilegeom,
                                  'west':west,'south':south,'east':east,'north':north}
       
    def _GetSystemDefRegions(self):
        '''
        '''
        #I need the layer for the region that is not in the regionsmodis table
        wherestatement = 'WHERE M.regionid IS NULL'
        recs = self.session._SelectAllDefRegions(wherestatement)

        for rec in recs:          
            compid = 'defaultregions_roi'
            folder = 'defaultregions'
            band = 'roi'
            regionid = rec[1]
            compD = {'compid': compid,'folder':folder, 'band':band}
            system = 'system'
            
            layerComp = self.session._SelectComp(system,compD)
            #
            system = 'system'
            comp = Composition(layerComp, system, self.process.system.srcdivision)
            #
            queryD = {}
            queryD['compid'] = comp.compid
            queryD['regionid'] = regionid.lower()

            paramL = ['source', 'product', 'suffix', 'acqdate', 'acqdatestr', 'doy', 'createdate', 'regionid']
            print ('finding mgrs for',rec)
            layerrec = self.session._SelectLayer(system,queryD,paramL)

            if layerrec == None:
                continue 
            source, product, suffix, acqdate, acqdatestr, doy, createdate, regionid = layerrec
            if acqdate == None:
                acqdate = False
            datumD = {'acqdatestr': acqdatestr,'acqdate':acqdate}
            #Set the locus         
            loc = regionid
            #Set the locuspath
            locusPath = regionid
            #Construct the locus dictionary
            locusD = {'locus':loc, 'locusPath':locusPath, 'path':locusPath}
            #Create the layer
            layer = VectorLayer(comp, locusD, datumD, self.process.srcpath)
            #Get the layer and the geom
            srcDS,srcLayer = mj_gis.ESRIOpenGetLayer(layer.FPN)
            for feature in srcLayer.layer:  
                geom = mj_gis.Geometry()
                #add the feature and extract the geom
                geom.GeomFromFeature(feature)
                if srcLayer.geomtype.lower() != 'polygon':
                    ERRORIGEN
                west, south, east, north = geom.shapelyGeom.bounds
                #Get the tiles inside this region
                tiles = self.session._SearchTilesFromWSEN(west, south, east, north)
                for tile in tiles:
                    #hvtile,htile,vtile,west,south,east,north,ullon,ullat,urlon,urlat,lrlon,lrlat,lllon,lllat, minx, miny, maxx, maxy = tile
                    mgrs,west,south,east,north,ullon,ullat,urlon,urlat,lrlon,lrlat,lllon,lllat, minx, miny, maxx, maxy = tile
                    llptL = ((ullon,ullat),(urlon,urlat),(lrlon,lrlat),(lllon,lllat))
                    tilegeom = mj_gis.Geometry()
                    tilegeom.PointsToPolygonGeom(llptL)
                    #Get the overlap
                    overlapGeom = tilegeom.ShapelyIntersection(self.senTileD[mgrs]['geom'])  
                    productoverlap = overlapGeom.area/tilegeom.shapelyGeom.area
                    if productoverlap >= 0:
                        query = {'system':'system', 'regionid':regionid,'regiontype':'default', 'overwrite':False, 'delete':False, 'mgrs':mgrs,'utmzone':self.senTileD[mgrs]['utmzone'], 'mgrsid':self.senTileD[mgrs]['mgrsid']}
                        self.session._InsertSentinelRegionTile(query)
  