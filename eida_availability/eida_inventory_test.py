#! /usr/bin/env python
#
# file eida_inventory_test.py
#      ======================
#
# K. Stammler, 17-Jun-2020

"""
Test 'get_stations' response, compare retrieved networks between routing
client and separate requests to EIDA servers directly.

Syntax:
    eida_inventory_test.py <level>
    
    <level>     request level, 'network', 'station' or 'channel'
"""

from __future__ import print_function
import os
import sys
import time
import datetime
from obspy import UTCDateTime, _get_version_string
from obspy.clients.fdsn import RoutingClient
from obspy.clients.fdsn.client import Client
#import pickle

eida_servers = (
    'NOA', 'http://eida.geo.uib.no', 'RESIF', 'BGR', 'ETH',
    'ODC', 'GFZ', 'INGV', 'NIEP', 'KOERI', 'LMU', 'ICGC',
)
wanted_channels = ( 'HHZ', 'BHZ', 'EHZ', 'SHZ' )
timeout = 240
legal_reqlevels = ('network','station','channel')

reference_networks = [
    'NL',   # ODC
    'GE',   # GFZ
    'FR',   # RESIF
    'CH',   # ETHZ
    'GR',   # BGR
    'BW',   # LMU
    'RO',   # NIEP
    'KO',   # KOERI
    'HL',   # NOA
    'NO',   # UIB/NORSAR
    'CA',   # ICGC
    'IV',   # INGV
]



#-------------------------------------------------------------------------------


class Logtext:
    logfile = os.path.join( os.environ['HOME'], 'Work', 'EidaAvailability',
        'eida_inventory_test.log' )
    def __init__( self ):
        if not os.path.exists(os.path.dirname(self.logfile)):
            os.makedirs( os.path.dirname(self.logfile) )
        self.fp = open( self.logfile, 'a' )
    def close( self ):
        self.fp.close()
    def write( self, text ):
        print( text )
        self.fp.write( "%s\n" % text ) 


#-------------------------------------------------------------------------------


def save_channel_inventory_for_eidaavailability( chaninv ):
    "Save channel level inventory for EidaAvailability program."
    import pickle
    eiainv = os.path.join( os.environ['HOME'], 'Work', 'EidaAvailability',
        'chanlist_cache.pickle' )
    if not os.path.exists(os.path.dirname(eiainv)):
        return 'no EidaAvailability working path'
    tmpfile = "/tmp/channelinv.pickle"
    fp = open( tmpfile, 'wb' )
    pickle.dump( chaninv, fp )
    fp.close()
    ret = os.system( "mv %s %s" % (tmpfile,eiainv) )
    if ret != 0:
        print( "*** channel inventory update failed ***" )
        return 'channel inventory save failed'
    else:
        return 'channel inventory saved for EidaAvailability.py'


#-------------------------------------------------------------------------------


def check_missing_reference_networks( rnets ):
    miss = []
    for net in reference_networks:
        if net not in rnets:
            miss.append( net )
    return sorted(miss)


#-------------------------------------------------------------------------------



if __name__ == '__main__':

    # Retrieve parameter from command line.
    if len(sys.argv) < 2:
        print( __doc__ )
        exit()
    reqlevel = sys.argv[1]
    if reqlevel not in legal_reqlevels:
        print( "Legal request levels are %s" % ','.join(legal_reqlevels) )
        exit()

    stamp = time.time()
    lt = Logtext()
    lt.write( "\neida_inventory_test.py started at %s MEST, level %s (obspy %s) timeout %d (timeout bugfix, no restricted)"
        % (datetime.datetime.now().strftime("%d-%b-%Y_%T"),reqlevel,
        _get_version_string(),timeout) )

    # Request interval is last year.
    t2 = UTCDateTime()
    t1 = t2 - 365*86400
    
    # 'get_stations' parameters.
    invpar = {
        'level'              : reqlevel,
        'channel'            : ','.join(wanted_channels),
        'starttime'          : t1,
        'endtime'            : t2,
        'includerestricted'  : False,
    }

    # Loop all EIDA servers for separate requests and store network list.
    snets = set( [] )
    for srv in eida_servers:
        lt.write( "    reading inventory from server %s" % srv )
        try:
            client = Client( srv, timeout=timeout )
            sinv = client.get_stations( **invpar )
        except Exception as e:
            lt.write( "        FAILED: %s" % repr(e) )
            continue
        addset = set( sinv.get_contents()['networks'] )
        snets = snets.union( addset )
    
    invpar['timeout'] = timeout
    
    # Use RoutingClient.
    lt.write( "    reading inventory from routing client" )
    try:
        roc = RoutingClient( "eida-routing" )
        rinv = roc.get_stations( **invpar )
    except Exception as e:
        lt.write( "        FAILED: %s" % repr(e) )
        lt.write( "\n==========================================================\n" )
        lt.close()
        exit()
    rnets = set( rinv.get_contents()['networks'] )
    missref = check_missing_reference_networks( rnets )
    
    runtime = time.time() - stamp

    # Write results.
    if missref:
        lt.write( "missing reference networks: %s" % ','.join(missref) )
    lt.write( "rnets (%d) %s" % (len(rnets),', '.join(sorted(rnets))) )
    lt.write( "snets (%d) %s" % (len(snets),', '.join(sorted(snets))) )
    lt.write( "rnets-snets %s" % ', '.join(sorted(rnets-snets)) )
    lt.write( "snets-rnets %s" % ', '.join(sorted(snets-rnets)) )
    lt.write( "runtime %3.1fs" % runtime )
    lt.write( "\n==========================================================\n" )
    
    #if len(rnets) > 110 and reqlevel == 'channel':
    #    msg = save_channel_inventory_for_eidaavailability( rinv )
    #    lt.write( "\n     eiasave: %s\n" % msg )

    lt.close()
