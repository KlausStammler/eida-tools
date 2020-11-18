#! /usr/bin/env python
#
# file MetadataWebRequest.py
#      ====================
#
# K. Stammler, 14-Dec-2017

"""
Generate and apply URL for FDSNws metadata request.

Syntax:
    MetadataWebRequest.py <net> [<station>] [<level>] [<format>]
    
    <net>       Net code like 'GR' (input converted to uppercase)
    <station>   Station name like 'GRA1' (input converted to uppercase)
                <station> must be part of <net>. If ommitted, lists all
                stations of network.
    <level>     Detail level. Possible values:
                'station':  one line per station (in text mode)
                'channel':  one line per channel (in text mode)
                'response': detailed response info (only xml format)
                (default: 'station' with empty <station> parameter,
                'channel' with <station> parameter)
    <format>    Output format. Either 'text' or 'xml', default 'text'.
                

Examples:
    MetadataWebRequest.py be
    MetadataWebRequest.py be rchb
    MetadataWebRequest.py be rchb response
    MetadataWebRequest.py be rchb channel xml
"""

from __future__ import print_function
import os
import sys


#-------------------------------------------------------------------------------


def metadata_request( net, station, oformat='text', level=None ):

    url = metadata_url( net, station, oformat='text', level=None )
    print( "requesting:", url )
    os.system( '%s "%s"' % (mwr_browser,url) )


#-------------------------------------------------------------------------------


def metadata_url( net, station, oformat='text', level=None, server=None ):

    if server:
        srv = server
    else:
        srv = get_server( net )
        if srv == None:
            print( "No server found for net '%s'" % net )
            exit()
    
    url = "%s/%s/" % (mwr_rootadr.get(srv,srv),mwr_station_ext)
    
    if level == None:
        if station == None:
            level = 'station'
        else:
            level = 'channel'
    if level == 'response':
        oformat = 'xml'
        
    url += "query?format=%s&level=%s&network=%s" % (oformat,level,net)
    
    if station != None:
        url += "&station=%s" % station
    return url


#-------------------------------------------------------------------------------


def parse_aliases( alias ):
    for srv in mwr_aliases.keys():
        if alias in mwr_aliases[srv]:
            return srv
    return alias


#-------------------------------------------------------------------------------


def get_server( net ):
    for srv in mwr_server.keys():
        if net in mwr_server[srv]:
            return srv
    #return None
    return 'IRIS'


#-------------------------------------------------------------------------------


def available_netcodes():
    nc = []
    for srv in mwr_server.keys():
        nc += mwr_server[srv]
    return sorted(nc)


#-------------------------------------------------------------------------------


mwr_rootadr = {
    'BGR' :   'http://eida.bgr.de/fdsnws',
    'ETH':    'http://eida.ethz.ch/fdsnws',
    'GFZ':    'http://geofon.gfz-potsdam.de/fdsnws',
    'ICGC':   'http://ws.icgc.cat/fdsnws',
    'IPGP':   'http://eida.ipgp.fr/fdsnws',
    'INGV':   'http://webservices.ingv.it/fdsnws',
    'IRIS':   'http://service.iris.edu/fdsnws',
    'KOERI':  'http://eida-service.koeri.boun.edu.tr/fdsnws',
    'LMU':    'http://erde.geophysik.uni-muenchen.de/fdsnws',
    'NIEP':   'http://eida-sc3.infp.ro/fdsnws',
    'NOA':    'http://eida.gein.noa.gr/fdsnws',
    'ODC':    'http://www.orfeus-eu.org/fdsnws',
    'RESIF':  'http://ws.resif.fr/fdsnws',
    'NCEDC':  'http://service.ncedc.org/fdsnws',
    'SCEDC':  'http://service.scedc.caltech.edu/fdsnws',
    'UIB'  :  'http://eida.geo.uib.no/fdsnws',
}

mwr_aliases = {
    'ODC' : ['ORFEUS'],
    'ETH' : ['ETHZ','SED'],
    'GFZ' : ['GEOFON'],
    'UIB' : ['NORSAR'],
}


mwr_server = {
    'BGR' : ['GR','RN','SX','TH'],
    'ETH' : ['CH','C4','S'],
    'INGV' : [
        'AC','BA','GU','IV','IX','MN','NI','OT','OX','RF','SI','ST','TV','XO'
    ],
    'GFZ' : [
        'AW','CK','CN','CX','CZ','DK','EE','EI','FN','GE','HE','HN','HT',
        'HU','IO','IS','JS','KC','KP','M1','NU','PL','PM','SJ','SK','TT','WM'
    ],
    'KOERI': ['IJ','KO','TL'],
    'LMU' : ['BW','Z3','ZJ'],
    'NIEP' : ['BS','MD','RO'],
    'NOA' : ['CQ','EG','HA','HC','HL','HP','ME'],
    'ODC' : [
        'AB','AI','BE','BN','CA','CR','DZ','EB','ES','GB','GO','HF','IB','II',
        'IP','LC','LX','NA','NL','NR','OE','SL','SS','TK','TU',
        'UP','VI','WC','YF'
    ],
    'RESIF' : ['CL','FR','MT','ND','RA','RD','G'],
    'NCEDC' : [
        '3B', '4B', '5B', '6B', 'AZ', 'BG', 'BP', 'CE', 'GM', 'GS',
        'LA', 'LB', 'NC', 'NN', 'NP', 'PB', 'PG', 'SF', 'TA', 'UL',
        'UO', 'US', 'UW', 'WR',
    ],
    'SCEDC' : ['BC', 'CI', 'MX', 'NC', 'SB', 'SN', 'ZY'],
    'IRIS' : [
        'BC', 'AE', 'AF', 'AG', 'AI', 'AK', 'AL', 'AM', 'AO', 'AR', 'AS', 'AT',
        'AU', 'AV', 'AY', 'BC', 'BF', 'BI', 'BL', 'BN', 'BX', 'C',  'C0', 'C1',
        'CA', 'CC', 'CD', 'IU',
    ],
    'UIB' : [ 'NO', 'NS' ],
}

mwr_station_ext = "station/1"

mwr_browser = "firefox"



if __name__ == '__main__':

    if len(sys.argv) < 2:
        print( __doc__ )
        print( "available netcodes: ", ','.join(available_netcodes()) )
        exit()
    
    net = sys.argv[1].upper()
    station = None
    level = None
    oformat = 'text'
    if len(sys.argv) > 2 and sys.argv[2]:
        station = sys.argv[2].upper()
    if len(sys.argv) > 3 and sys.argv[3]:
        if sys.argv[3] in ('station','channel','response'):
            level = sys.argv[3]
    if len(sys.argv) > 4 and sys.argv[4]:
        if sys.argv[4] in ('text','xml'):
            oformat=sys.argv[4]
    
    metadata_request( net, station, oformat=oformat, level=level )
