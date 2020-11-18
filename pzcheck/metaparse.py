#! /usr/bin/env python
#
# file metaparse.py
#      ============
#
# K. Stammler, 7-Dec-2019

"""
Parse station metadata XML read from file or from FDSN server.

Syntax:
    metaparse.py <cmd> <xmlfile>
    metaparse.py <cmd> <netsta> [<server>]

    <cmd>: 'c' or 'compact'   : compact print
           'e' or 'extended'  : extended output
           'p' or 'plaindump' : plain dump
           'x' or 'xml'       : XML output (copy of file or server response)
           'n' or 'normcheck' : check poles&zeros in each stage for:
                                - normalization at given frequency
                                - complex conjugate pairs
                                - negative real part
    <netsta>  : network and station, e.g. 'GR.BFO', the station part accepts
                wildcards (use with parantheses)
    <server>  : server to access, default is default server for network

Examples:
    metaparse.py c gr.gra1
    metaparse.py n 'gr.*'
"""

from __future__ import print_function
import os
import sys
if sys.version_info.major == 3:
    from urllib.request import urlopen
    from io import StringIO
else:
    from urllib2 import urlopen
    from StringIO import StringIO
from lxml import etree
import numpy as np
from scipy import signal
import MetadataWebRequest


normdevtolerance = 2.0  # in percent

#-------------------------------------------------------------------------------


def csort( carr ):
    "Returns sorted complex array."
    return np.array(sorted(np.array(carr)))


#-------------------------------------------------------------------------------


pzlib = {

    'STS2-GRSN-simple' : [
        csort([(-0.0367429+0.036754j), (-0.0367429-0.036754j)]),
        csort([0j, 0j]),
    ],
    'STS1-GRF' : [
        csort([(-0.222111+0.222178j), (-0.222111-0.222178j),
            (-31.416+0j), (-19.572+24.574j), (-19.572-24.574j),
            (-7.006+30.625j), (-7.006-30.625j), (-28.306+13.629j),
            (-28.306-13.629j)]),
        csort([0j, 0j]),
    ],
    'STS2-generic' : [
        csort([(-0.037004+0.037016j), (-0.037004-0.037016j), (-251.33+0j),
            (-131.04-467.29j), (-131.04+467.29j)]),
        csort([0j, 0j]),
    ],
    'STS2-gen3' : [
        csort([(-13300+0j), (-10530+10050j), (-10530-10050j), (-520.3+0j),
            (-374.8+0j), (-97.34+400.7j), (-97.34-400.7j), (-15.64+0j),
            (-0.037+0.037j), (-0.037-0.037j), (-255.1+0j)]),
        csort([0j, 0j, (-463.1+430.5j), (-463.1-430.5j), (-176.6+0j),
            (-15.15+0j)]),
    ],
    'STS1-GRF-simple' : [
        csort([(-0.222111+0.222178j), (-0.222111-0.222178j)]),
        csort([0j, 0j]),
    ],
    'CMG-3T' : [
        csort([(-0.03701+0.03701j), (-0.03701-0.03701j), (-1131+0j),
            (-1005+0j), (-502.7+0j)]),
        csort([0j, 0j]),
    ],
    'CMG-3TB-360' : [
        csort([(-0.01234+0.01234j), (-0.01234-0.01234j), (-188.828+195.54j),
            (-188.828-195.54j), (-259.222+719.645j), (-259.222-719.645j)]),
        csort([0j, 0j]),
    ],
    'CMG-3ESPC_60s' : [
        csort([(-0.074016+0.074016j), (-0.074016-0.074016j), (-502.65+0j),
            (-1005+0j), (-1131+0j)]),
        csort([0j, 0j]),
    ],
    'CMG-3ESP_120s' : [
        csort([(-0.037008+0.037008j), (-0.037008-0.037008j), (-502.65+0j),
            (-1005+0j), (-1131+0j)]),
        csort([0j, 0j]),
    ],
    'Meridian-PH/120' : [
        csort([(-0.036614+0.037059j), (-0.036614-0.037059j), (-32.55+0j),
            (-142+0j), (-364+404j), (-364-404j), (-1260+0j), (-4900+5200j),
            (-4900-5200j), (-7100+1700j), (-7100-1700j)]),
        csort([0j, 0j, (-31.63+0j), (-160+0j), (-350+0j), (-3177+0j)]),
    ],
    'Trillium-120' : [
        csort([(-0.03859+0.03649j), (-0.03859-0.03649j), (-190+0j),
            (-158+193j), (-158-193j), (-639+1418j), (-639-1418j)]),
        csort([0j, 0j, (-106+0j), (-158+0j)]),
    ],
    'TrilliumCompact-20' : [
        csort([(-0.2214+0.2221j), (-0.2214-0.2221j), (-343+0j), (-370+467j),
            (-370-467j), (-836+1522j), (-836-1522j), (-4900+4700j),
            (-4900-4700j), (-6900+0j), (-15000+0j)]),
        csort([0j, 0j, (-392+0j), (-1960+0j), (-1490+1740j), (-1490-1740j)]),
    ],
    'KS-36000-I' : [
        csort([(-89.85+0j), (-18.43+18.91j), (-18.43-18.91j),
            (-0.01234+0.01234j), (-0.01234-0.01234j), (-0.004219+0j)]),
        csort([0j, 0j, 0j]),
    ],
    'Lennartz_LE-3Dlite' : [
        csort([(-4.44+4.44j), (-4.44-4.44j), (-1.083+0j)]),
        csort([0j, 0j, 0j]),
    ],
    'Lennartz_LE-3D5s' : [
        csort([(-0.888+0.888j), (-0.888-0.888j), (-0.29+0j)]),
        csort([0j, 0j, 0j]),
    ],
    'GS13' : [
        csort([(-4.443+4.443j), (-4.443-4.443j)]),
        csort([0j, 0j])
    ],
    'Vegik-SP' : [
        csort([(-13.6006+33.2863j), (-13.6006-33.2863j), (-21.175+24.2739j),
            (-21.175-24.2739j), (-25.5848+15.9695j), (-25.5848-15.9695j),
            (-27.9812+7.9335j), (-27.9812-7.9335j), (-28.7466+0j),
            (-0.037004+0.037016j), (-0.037004-0.037016j)]),
        csort([0j, 0j]),
    ],
}


#-------------------------------------------------------------------------------


class MetaXML:

    def __init__( self, xmltext ):
        self.xmltext = xmltext
        self.header = {}
        self.stainfo = {}
        self.xmlroot = self.check_and_get_xml_root()
        self.parse_header()
        self.parse_tree()
    
    def check_and_get_xml_root( self ):
        if self.xmltext.strip() == '':
            print( "FDSN/EIDA service returns empty XML text" )
            exit()
        try:
            xmlroot = etree.fromstring(self.xmltext)
        except:
            print( self.xmltext )
            print( "Parser error in XML code" )
            exit()
        return xmlroot
    
    def printmeta( self, cmd ):
        if cmd in ('c','compact'):
            self.print_compact()
        elif cmd in ('e','extended'):
            self.dump()
        elif cmd in ('p','plaindump'):
            self.plaindump()
        elif cmd in ('x','xml'):
            print( self.xmltext )
        elif cmd in ('n','checknorm'):
            self.check_pz_norm()
        else:
            print( "unknown cmd '%s'" % cmd )
    
    def parse_header( self ):
        for rootelem in self.xmlroot.getchildren():
            if tagstrip(rootelem.tag) != "Network":
                self.header[tagstrip(rootelem.tag)] = rootelem.text
    
    def get_station_nodes( self ):
        netnodes = []
        for rootelem in self.xmlroot.getchildren():
            if tagstrip(rootelem.tag) != "Network":
                continue
            netnodes.append( rootelem )
        if netnodes == []:
            return []
        stanodes = []
        for netnode in netnodes:
            netcode = netnode.get('code')
            for netelem in netnode.getchildren():
                if tagstrip(netelem.tag) == "Station":
                    stanodes.append( (netcode,netelem) )
        return stanodes
    
    def parse_tree( self ):
        for netcode,stanode in self.get_station_nodes():
            self.parse_stanode( netcode, stanode )
    
    def parse_stanode( self, netcode, stanode ):
        netsta = "%s.%s" % (netcode,stanode.get('code'))
        if not netsta in self.stainfo:
            self.stainfo[netsta] = []
        stanodeinfo = {
            'startDate' : stanode.get('startDate'),
            'endDate' : stanode.get('endDate'),
        }
        channels = {}
        chancnt = {}
        for staelem in stanode.getchildren():
            chancode = staelem.get('code')
            if chancode in chancnt.keys():
                chancnt[chancode] += 1
            else:
                chancnt[chancode] = 0
            chanid = "%s-%d" % (chancode,chancnt[chancode])
            for elemname in ('Latitude','Longitude','Elevation','Depth'):
                if tagstrip(staelem.tag) == elemname:
                    stanodeinfo[elemname] = float( staelem.text )
            if tagstrip(staelem.tag) == 'Channel':
                channels[chanid] = self.parse_channel( staelem )
        stanodeinfo['channels'] = channels
        self.stainfo[netsta].append( stanodeinfo )
    
    def parse_channel( self, staelem ):
        chaninfo = {
            'startDate' : staelem.get('startDate'),
            'endDate' : staelem.get('endDate'),
        }
        for chanelem in staelem.getchildren():
            for elemname in ('Latitude','Longitude','Elevation','Depth',
                'Azimuth','Dip'):
                if tagstrip(chanelem.tag) == elemname:
                    chaninfo[elemname] = float( chanelem.text )
            if tagstrip(chanelem.tag) == 'SampleRateRatio':
                chaninfo['SampleRate'] = self.parse_sample_rate( chanelem )
            for elemname in ('Sensor','DataLogger'):
                if tagstrip(chanelem.tag) == elemname:
                    chaninfo[elemname] = self.parse_device( chanelem )
            if tagstrip(chanelem.tag) == 'Response':
                chaninfo['Response'] = self.parse_response( chanelem )
        return chaninfo
    
    def parse_response( self, respelem ):
        resp = {}
        stages = {}
        for elem in respelem.getchildren():
            if tagstrip(elem.tag) == 'InstrumentSensitivity':
                resp['TotalSens'], resp['SensFreq'], resp['SensUnits'] = \
                    self.parse_instrument_sensitivity( elem )
            if tagstrip(elem.tag) == 'Stage':
                stages[elem.get('number')] = self.parse_stage( elem )
        resp['stages'] = []
        for stagenum in sorted(stages.keys()):
            resp['stages'].append( stages[stagenum] )
        return resp
    
    def parse_stage( self, stagelem ):
        stageinfo = {}
        for elem in stagelem.getchildren():
            if tagstrip(elem.tag) == 'StageGain':
                stageinfo['StageGainValue'], stageinfo['StageGainFrq'] = \
                    self.parse_stage_gain( elem )
            elif tagstrip(elem.tag) == 'PolesZeros':
                stageinfo['PolesZeros'] = self.parse_poles_zeros( elem )
            elif tagstrip(elem.tag) == 'Decimation':
                stageinfo['Decimation'] = self.parse_decimation( elem )
            elif tagstrip(elem.tag) == 'FIR':
                stageinfo['FIR'], stageinfo['FIRSymmetry'] \
                    = self.parse_fir( elem )
        return stageinfo
    
    def parse_decimation( self, decelem ):
        decinfo = {}
        for elem in decelem.getchildren():
            for tagname in ('InputSampleRate','Factor','Offset','Delay',
                'Correction'):
                if tagstrip(elem.tag) == tagname:
                    decinfo[tagname] = float( elem.text )
        return decinfo
    
    def parse_fir( self, firelem ):
        symmetry = 'unknown'
        fircoef = []
        for elem in firelem.getchildren():
            if tagstrip(elem.tag) == 'Symmetry':
                symmetry = elem.text
            elif tagstrip(elem.tag) == 'NumeratorCoefficient':
                fircoef.append( float(elem.text) )
        return (fircoef,symmetry)
    
    def parse_poles_zeros( self, pzelem ):
        pzinfo = {}
        pzinfo['units'] = 'unknown'
        pzinfo['poles'] = []
        pzinfo['zeros'] = []
        inputunit = None
        outputunit = None
        for elem in pzelem.getchildren():
            if tagstrip(elem.tag) == 'InputUnits':
                child = elem.getchildren()[0]
                inputunit = child.text
            elif tagstrip(elem.tag) == 'OutputUnits':
                child = elem.getchildren()[0]
                outputunit = child.text
            elif tagstrip(elem.tag) == 'PzTransferFunctionType':
                pzinfo['pztype'] = elem.text.lower()
            elif tagstrip(elem.tag) == 'NormalizationFactor':
                pzinfo['norm'] = float( elem.text )
            elif tagstrip(elem.tag) == 'NormalizationFrequency':
                pzinfo['normfrq'] = float( elem.text )
            elif tagstrip(elem.tag) == 'Zero':
                pzinfo['zeros'].append( self.parse_complex_number(elem) )
            elif tagstrip(elem.tag) == 'Pole':
                pzinfo['poles'].append( self.parse_complex_number(elem) )
        if inputunit and outputunit:
            pzinfo['units'] = "%s -> %s" % (inputunit.lower(),outputunit.lower())
        return pzinfo
    
    def parse_complex_number( self, celem ):
        r = 0.
        i = 0.
        for elem in celem:
            if tagstrip(elem.tag) == 'Real':
                r = float( elem.text )
            elif tagstrip(elem.tag) == 'Imaginary':
                i = float( elem.text )
        return complex( r, i )
    
    def parse_stage_gain( self, sgelem ):
        value = 1.0
        frq = 1.0
        for elem in sgelem.getchildren():
            if tagstrip(elem.tag) == 'Value':
                value = float( elem.text )
            elif tagstrip(elem.tag) == 'Frequency':
                frq = float( elem.text )
        return (value,frq)

    def parse_instrument_sensitivity( self, senselem ):
        totsens = 1.
        sensfrq = 1.
        sensunits = 'unknown'
        inputunit = None
        outputunit = None
        for elem in senselem.getchildren():
            if tagstrip(elem.tag) == 'Value':
                totsens = float( elem.text )
            elif tagstrip(elem.tag) == 'Frequency':
                sensfrq = float( elem.text )
            elif tagstrip(elem.tag) == 'InputUnits':
                child = elem.getchildren()[0]
                inputunit = child.text
            elif tagstrip(elem.tag) == 'OutputUnits':
                child = elem.getchildren()[0]
                outputunit = child.text
        if inputunit and outputunit:
            sensunits = "%s -> %s" % (inputunit.lower(),outputunit.lower())
        return (totsens,sensfrq,sensunits)

    def parse_sample_rate( self, srelem ):
        numsamples = 1
        numseconds = 1
        for elem in srelem.getchildren():
            if tagstrip(elem.tag) == 'NumberSamples':
                numsamples = int( elem.text )
            elif tagstrip(elem.tag) == 'NumberSeconds':
                numseconds = int( elem.text )
        return float(numsamples)/float(numseconds)
    
    def parse_device( self, develem ):
        devnames = []
        for elem in develem.getchildren():
            for elemname in ('Type','Description','Manufacturer','Model'):
                if tagstrip(elem.tag) == elemname:
                    devnames.append( elem.text )
        return devnames
    
    def print_header( self ):
        for k in sorted(self.header.keys()):
            print( k, self.header[k] )
    
    def plaindump( self ):
        print( self.stainfo )
        
    def dump( self ):
        self.print_header()
        #elemhide = ['Azimuth','Elevation','Dip','Depth']
        elemhide = []
        chanhide = ['EX1','EX2','EX3','GAN','GEL','GLA','GLO','GNS','GPL',
            'GST','LCE','LCQ','VCO','VDT','VEC','VEI','VM1','VM2','VM3','VPB']
        for sta in self.stainfo.keys():
            print( sta )
            for staelem in self.stainfo[sta]:
                print()
                for kc in sorted(staelem.keys()):
                    if kc in elemhide: continue
                    if kc != 'channels':
                        print( '  ', kc, staelem[kc] )
                for chan in sorted(staelem['channels'].keys()):
                    if chan in chanhide: continue
                    print( '    ', chan )
                    for kd in sorted(staelem['channels'][chan].keys()):
                        if kd in elemhide: continue
                        if kd != 'Response':
                            print( '      ', kd, staelem['channels'][chan][kd] )
                    print( '      ', 'Response' )
                    for ks in sorted(staelem['channels'][chan]['Response'].keys()):
                        if ks in elemhide: continue
                        if ks == 'TotalSens':
                            totsens = staelem['channels'][chan]['Response'][ks]
                            print( '        ', ks, totsens, 'or', 1e9/totsens )
                        elif ks != 'stages':
                            print( '        ', ks, staelem['channels'][chan]['Response'][ks] )
                    print( '         stages' )
                    for snum,stage in enumerate(staelem['channels'][chan]['Response']['stages']):
                        print( '          ', snum+1, ':' )
                        for kst in sorted(stage.keys()):
                            if kst in elemhide: continue
                            if kst != 'PolesZeros':
                                print( '            ', kst, stage[kst] )
                        if 'PolesZeros' in stage.keys():
                            print( '             PolesZeros:' )
                            for kpz in sorted(stage['PolesZeros'].keys()):
                                if kpz in elemhide: continue
                                print( '              ', kpz, stage['PolesZeros'][kpz] )

    def print_compact( self ):
        chanhide = ['EX1','EX2','EX3','GAN','GEL','GLA','GLO','GNS','GPL',
            'GST','LCE','LCQ','VCO','VDT','VEC','VEI','VM1','VM2','VM3','VPB']
        for sta in self.stainfo.keys():
            staline = "%s" % sta
            stalat = None
            stalon = None
            staelev = None
            stadepth = 0.
            if 'Latitude' in self.stainfo[sta][0].keys()\
                and 'Longitude' in self.stainfo[sta][0].keys():
                stalat = self.stainfo[sta][0]['Latitude']
                stalon = self.stainfo[sta][0]['Longitude']
                staline += "    (%8.5f,%8.5f)" % (stalat,stalon)
            if 'Depth' in self.stainfo[sta][0].keys():
                stadepth = self.stainfo[sta][0]['Depth']
            if 'Elevation' in self.stainfo[sta][0].keys():
                staelev = self.stainfo[sta][0]['Elevation']
                staline += " at %5.1fm" % (staelev-stadepth)
            print( staline )
            for staelem in self.stainfo[sta]:
                validinterv = self.pcmp_get_validinterv( staelem )
                print( "  %s" % validinterv )
                altstalat = "mainlat"
                altstalon = "mainlon"
                altstaelev = "mainelev"
                altdepth = 0.
                if 'Latitude' in staelem.keys() \
                    and staelem['Latitude'] != stalat:
                    altstalat = "%g" % staelem['Latitude']
                if 'Longitude' in staelem.keys() \
                    and staelem['Longitude'] != stalon:
                    altstalon = "%g" % staelem['Longitude']
                if 'Depth' in staelem.keys():
                    altdepth = staelem['Depth']
                if 'Elevation' in staelem.keys() \
                    and staelem['Elevation'] != staelev:
                    altstaelev = "%g" % (staelem['Elevation']-altdepth)
                if altstalat != "mainlat" or altstalon != "mainlon":
                    print( "  **alt lat&lon  %s,%s" % (altstalat,altstalon) )
                if altstaelev != "mainelev":
                    print( "  **alt elev %s" % altstaelev )
                clprint = LineStore()
                for chan in sorted(staelem['channels'].keys()):
                    if chan in chanhide: continue
                    chanlines = []
                    chandict = staelem['channels'][chan]
                    chanheader = "    %s   %s  (%5.1fHz)" % (chan,
                        self.pcmp_get_validinterv(chandict),
                        chandict.get('SampleRate',0.0))
                    altstalat = "mainlat"
                    altstalon = "mainlon"
                    altstaelev = "mainelev"
                    altstadepth = 0.
                    if 'Latitude' in chandict.keys() \
                        and chandict['Latitude'] != stalat:
                        altstalat = "%g" % chandict['Latitude']
                    if 'Longitude' in chandict.keys() \
                        and chandict['Longitude'] != stalon:
                        altstalon = "%g" % chandict['Longitude']
                    if 'Depth' in chandict.keys():
                        altstadepth = chandict['Depth']
                    if 'Elevation' in chandict.keys() \
                        and chandict['Elevation'] != staelev:
                        alststaelev = "%g" % (chandict['Elevation']-altstadepth)
                    if altstalat != "mainlat" or altstalon != "mainlon":
                        chanlines.append( "      **alt lat&lon  %s,%s" % (
                            altstalat,altstalon) )
                    if altstaelev != "mainelev":
                        chanlines.append( "      **alt elev %s" % altstaelev )
                    elif 'Depth' in chandict.keys():
                        if altstaelev == "mainelev":
                            if staelev is None:
                                height = "unknown"
                            else:
                                height = "%5.1f" % (staelev - chandict['Depth'])
                        else:
                            height = "%5.1f" % (altstaelev - chandict['Depth'])
                        chanlines.append( "      abs height: %s" % height )
                    chanlines.append( "      Sensor: %s  Logger: %s" % (
                        self.pcmp_print_device(chandict.get('Sensor',None)),
                        self.pcmp_print_device(chandict.get('DataLogger',None),
                            dev='Logger')) )
                    chanlines += self.pcmp_print_response(
                        chandict.get('Response',None) )
                    chantext = clprint.print_lines( chan, chanlines )
                    if chantext.count('\n') == 0 and len(chantext) < 15:
                        print( chanheader + chantext )
                    else:
                        print( chanheader )
                        print( chantext )
    
    def pcmp_get_validinterv( self, elem ):
        sdate = 'now       '
        edate = 'now       '
        if 'startDate' in elem.keys() and elem['startDate']:
            sdate = elem['startDate'].split('T')[0]
        if 'endDate' in elem.keys() and elem['endDate']:
            edate = elem['endDate'].split('T')[0]
        return "%s to %s" % (sdate,edate)
    
    def pcmp_print_device( self, wordarr, dev='Sensor' ):
        if not wordarr:
            return 'unknown'
        if len(wordarr) == 1:
            return wordarr[0]
        if dev == 'Sensor':
            keyw = ('Streckeisen','STS-2','TrilliumCompact','Hyperion',
                'Lennartz','LE-3D','GS13','Guralp','ESP-3V','CMG-3TB','Mark',
                'Nanometrics','Meridian','Trillium-120','Geotech','KS-36000')
        else:
            keyw = ('Quanterra','Q330','EarthData','Earth Data','Reftek-130',
                'Reftek130','Nanometrics','HRD24','Guralp')
        wordscore = {}
        maxscore = 0
        secmaxscore = 0
        for word in wordarr:
            wordscore[word] = 0
            for w in keyw:
                if w.lower() in word.lower():
                    wordscore[word] += 1
            if wordscore[word] > maxscore:
                secmaxscore = maxscore
                maxscore = wordscore[word]
            elif wordscore[word] > secmaxscore:
                secmaxscore = wordscore[word]
        if maxscore == 0:
            return repr(wordarr)
        selwords = []
        for word in wordscore.keys():
            if wordscore[word] == maxscore:
                selwords.append( word )
        if len(selwords) < 2 and secmaxscore > 0:
            for word in wordscore.keys():
                if wordscore[word] == secmaxscore:
                    selwords.append( word )
        if len(selwords) < 2:
            return repr(wordarr)
        while len(selwords) > 2:
            maxlen = 0
            for w in selwords:
                if len(w) > maxlen:
                    maxlen = len(w)
            for w in selwords:
                if len(w) == maxlen:
                    selwords.remove( w )
                    break
        if len(selwords[0]) > len(selwords[1]):
            word1 = selwords[0]
            word2 = selwords[1]
        else:
            word1 = selwords[1]
            word2 = selwords[0]
        if word2.lower() in word1.lower():
            return word1
        else:
            return "%s %s" % (word1,word2)
    
    def pcmp_print_response( self, resp ):
        if not resp:
            return ["no reponse specified"]
        rlines = []
        if 'SensUnits' in resp.keys() and resp['SensUnits'] != "m/s -> counts":
            rlines.append( "SensUnits: %s" % resp['SensUnits'] )
        if 'TotalSens' in resp.keys():
            totsensstr = "%g nm/s/count" % (1e9/resp['TotalSens'])
        else:
            totsensstr = "unknown"
        if 'SensFreq' in resp.keys():
            sensfreqstr = "%5.2f Hz" % resp['SensFreq']
        else:
            sensfreqstr = "unknown"
        rlines.append( "sens %s @ %s" % (totsensstr,sensfreqstr) )
        if 'stages' in resp.keys():
            for stage in resp['stages']:
                # Ignore decimation stages (FIR)
                if 'Decimation' in stage.keys():
                    continue
                if 'StageGainValue' in resp.keys():
                    gainstr = "gain: %g" % resp['StageGainValue']
                    if 'StageGainFrq' in resp.keys():
                        gainstr += " @ %5.2f Hz" % resp['StageGainFrq']
                    rlines.append( gainstr )
                rlines += self.pcmp_pz_lines( stage['PolesZeros'] )
        return ['      %s' % ln for ln in rlines]
    
    def pcmp_pz_lines( self, pzinfo ):
        pzlines = []
        if 'norm' in pzinfo.keys():
            normline = "norm: %g" % pzinfo['norm']
            if 'normfrq' in pzinfo.keys():
                normline += " @ %5.2f Hz" % pzinfo['normfrq']
            pzlines.append( normline )
        if 'units' in pzinfo.keys():
            if pzinfo['units'] != "m/s -> v":
                pzlines.append( "units: %s" % pzinfo['units'] )
        if 'poles' in pzinfo.keys():
            poles = pzinfo['poles']
        else:
            poles = []
        if 'zeros' in pzinfo.keys():
            zeros = pzinfo['zeros']
        else:
            zeros = []
        pzname = self.check_pzlib( poles, zeros )
        if pzname:
            pzlines.append( "pzlib: %s" % pzname )
        else:
            pzlines.append( "poles: %s" % repr(poles) )
            pzlines.append( "zeros: %s" % repr(zeros) )
        return pzlines
    
    def check_pzlib( self, poles, zeros, maxrms=1.e-5 ):
        for pzname in pzlib.keys():
            libpoles,libzeros = pzlib[pzname]
            if len(libpoles) != len(poles) or len(libzeros) != len(zeros):
                continue
            if len(libpoles) == 0:
                poldif = 0.
            else:
                poldif = abs( libpoles - csort(poles) )
            prms = np.sqrt( (poldif*poldif.conj()).sum() )
            if len(libzeros) == 0:
                zerodif = 0.
            else:
                zerodif = abs( libzeros - csort(zeros) )
            zrms = np.sqrt( (zerodif*zerodif.conj()).sum() )
            rms = prms + zrms
            if rms < maxrms:
                return pzname
        return None
    
    def check_pz_norm( self ):
        for netsta in self.stainfo.keys():
            for elem in self.stainfo[netsta]:
                for idict in elem.keys():
                    if idict != 'channels':
                        continue
                    for chan in elem[idict]:
                        for celem in elem[idict][chan].keys():
                            if celem != 'Response':
                                continue
                            for relem in elem[idict][chan][celem].keys():
                                if relem != 'stages':
                                    continue
                                for stage in elem[idict][chan][celem][relem]:
                                    for selem in stage.keys():
                                        if selem != 'PolesZeros':
                                            continue
                                        self.check_stage_norm( netsta, chan,
                                            stage[selem] )

    def check_stage_norm( self, netsta, chan, pzdict ):
        status = ""
        infostr = "%s %s" % (netsta,chan)
        cjstr = self.check_complex_conjugates_and_sign( pzdict['zeros'] )
        if cjstr:
            status += "    zeros: %s\n" % cjstr
        cjstr = self.check_complex_conjugates_and_sign( pzdict['poles'], True )
        if cjstr:
            status += "    poles: %s\n" % cjstr
        needkeys = ('normfrq','norm','zeros','poles')
        misskeys= []
        for key in needkeys:
            if key not in pzdict.keys():
                misskeys.append( key )
        if misskeys != []:
            print( "%s: missing stage info: %s" % (infostr,','.join(misskeys)) )
            return
        pz = signal.ZerosPolesGain( pzdict['zeros'], pzdict['poles'],
             pzdict['norm'] )
        frq, amp = pz.freqresp( [2.*np.pi*pzdict['normfrq']] )
        ampres = amp[0]
        onedev = (abs(ampres)-1.)*100.0
        if abs(onedev) > normdevtolerance:
            status += "    normdev %3.1f%%" % onedev
        if status:
            print( "%s\n%s" % (infostr,status.rstrip('\n')) )
            print( "    pzdump:", pzdict )
        else:
            print( "%s OK" % infostr )
    
    def check_complex_conjugates_and_sign( self, zlist, ispoles=False ):
        "Returns error string or empty string if ok"
        czlist = zlist[:]
        while len(czlist) > 0:
            ztest = czlist[0]
            if ispoles and ztest.real > 0.:
                return "non-negative real part in %g+%gj" % (ztest.real,ztest.imag)
            if self.isreal(ztest):
                czlist = czlist[1:]
                continue
            newczlist = []
            for z in czlist[1:]:
                if ztest is not None and self.isreal(ztest*z):
                    ztest = None
                else:
                    newczlist.append( z )
            if ztest is not None:
                return "not paired: %g+%gj" % (ztest.real,ztest.imag)
            czlist = newczlist[:]
        return ""
    
    def isreal( self, z ):
        if abs(z.real) < 1.e-8:
            return (abs(z.imag) <= abs(z.real))
        zrat = abs(z.imag/z.real)
        return (zrat < 1.e-8)

#-------------------------------------------------------------------------------


class LineStore:

    def __init__( self ):
        self.store = {}
    
    def print_lines( self, name, lines ):
        if name in self.store.keys():
            return '\n'.join(lines)
        res = self.check_lines( lines )
        if res is None:
            self.store[name] = lines
            return '\n'.join(lines)
        else:
            return "      -> %s" % res
    
    def check_lines( self, lines ):
        for k in self.store.keys():
            if len(self.store[k]) != len(lines):
                continue
            for line1,line2 in zip(self.store[k],lines):
                if line1 != line2:
                    break
            else:
                return k
        return None


#-------------------------------------------------------------------------------


def tagstrip( tag ):
    if '}' in tag:
        return tag.split('}')[-1]
    else:
        return tag


#-------------------------------------------------------------------------------


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print( __doc__ )
        exit()
    
    cmd = sys.argv[1]

    server = None
    if len(sys.argv) > 3 and sys.argv[3]:
        server = sys.argv[3]
    
    xmlfile = sys.argv[2]
    if os.path.exists(xmlfile):
        with open(xmlfile) as fobj:
            xmltext = fobj.read()
        MetaXML( xmltext ).printmeta(cmd)
    else:
        try:
            net, sta = xmlfile.upper().split('.')
        except:
            print( __doc__ )
            exit()
        url = MetadataWebRequest.metadata_url( net, sta, oformat='xml',
            level='response', server=server )
        MetaXML( urlopen(url).read() ).printmeta(cmd)
    
