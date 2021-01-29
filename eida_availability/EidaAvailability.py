#! /usr/bin/env python
#
# file EidaAvailability.py
#      ===================
#
# K. Stammler, 17-Jun-2020
#
# For creation of reports, needs 'pandoc' and 'convert' (imagemagick).
# This also requires implementation specific code, see variables eia_spec_...
#
# Create report with:
# import EidaAvailability
# ar = EidaAvailability.AvailabilityReport()
# ar.make_md_report( 'avreport.md' )
# pdfreport = ar.make_pdf_report( 'avreport.md' )

"""
Availability test of EIDA stations using Python obspy library.

- Conducts random waveform requests to single channels of EIDA stations.
- One request per minute.
- Requested time span randomly selected from last year, span length between
  60 and 600 s.
- Station randomly selected from the subset of unrestricted European EIDA
  stations offering at least one out of channels `HHZ`, `BHZ`, `EHZ` or `SHZ`.
- Request full station metadata from selected station and choose channel
  randomly, restricted to channels `HH?`, `BH?`, `EH?` and `SH?`.
- On successful request apply a restitution to the waveform data.
- Evaluate and store result of request in a file database.
- Plot and statistically analyze content of file database.

The code does not use the waveform catalog, therefore empty waveform returns
are due to data gaps or due to problems in data access and delivery.
"""


from __future__ import print_function
import os
import sys
import time
import signal
import datetime
import pickle
import numpy as np
from obspy.clients.fdsn import RoutingClient
from obspy import UTCDateTime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, projection_params

eia_global_timespan_days = 365  # Waveform selections within this day span.
eia_timeout = 60                # Timeout for retrieving station metadata.
eia_datapath = os.path.join( os.environ['HOME'], 'Work', 'EidaAvailability' )
eia_min_num_networks = 80       # With uncluderestricted=False.
eia_reqstats_timespan_days = 92 # Request statistics over 3 months back.

# Implementation specific settings for creating reports:
eia_spec_default_cssfile = "/home/stammler/svn/SzgrfDc/texts/styles/ks.css"
eia_spec_invtest_host = "detector@amp2f"  # Here the inventory test is run.
eia_spec_invtest_rempath = "/home/detector/Work/EidaAvailability"

# Status codes
STATUS_OK = 0
STATUS_NOSERV = 1
STATUS_METAFAIL = 2
STATUS_NODATA = 3
STATUS_FRAGMENT = 4
STATUS_INCOMPLETE = 5
STATUS_RESTFAIL = 6

# Names of status codes
error_names = {
    STATUS_OK         : 'OK',
    STATUS_NOSERV     : 'NOSERV',
    STATUS_METAFAIL   : 'METAFAIL',
    STATUS_NODATA     : 'NODATA',
    STATUS_FRAGMENT   : 'FRAGMENT',
    STATUS_INCOMPLETE : 'INCOMPL',
    STATUS_RESTFAIL   : 'RESTFAIL',
}

# Exclude temporary and not european networks.
exclude_networks = [
    '1N', '1T', '3C', '4H', '5M', '7A', '8C', '9C', '9H', 'XK',
    'XN', 'XT', 'XW', 'YW', 'YZ', 'Z3', 'ZF', 'ZJ', 'ZM', 'ZS',
]
not_european = [
    'AI', 'AW', 'CK', 'CN', 'CX', 'GL', 'IO', 'IQ', 'KC', 'KP', 'MQ',
    'NA', 'ND', 'NU', 'PF', 'WC', 'WI'
]
exclude_networks += not_european

# Reduce probability for selection on large dense networks.
large_networks = {
    'NL' : 0.5,
}

# To check whether all servers responded to channel level full
# EIDA metadata request.
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


class EidaAvailability:

    wanted_channels = ( 'HHZ', 'BHZ', 'EHZ', 'SHZ' )
    global_span = ( UTCDateTime()-86400*eia_global_timespan_days,
        UTCDateTime() )
    maxcacheage = 5*86400  # Renew inventory file every 5 days
    minreqlen = 60         # Waveform request window, minimum length ...
    maxreqlen = 600        # ... and maximum length
    
    def __init__( self ):
        self.slist_cache = os.path.join( eia_datapath, 'chanlist_cache.pickle' )
        self.roc = RoutingClient( "eida-routing" )
        self.meta_time = None
        self.wave_time = None
        self.status = None
        self.requestpar = None
        self.trymgr = RetryManager( 'eidainventory', 4800 )
        self._check_path()
    
    def _check_path( self ):
        if not os.path.exists(eia_datapath):
            os.makedirs( eia_datapath )

    def _get_inventory_from_service( self ):
        "Retrieve station inventory from EIDA routing client."
        try:
            slist = self.roc.get_stations( level='channel',
                channel=','.join(self.wanted_channels),
                starttime=self.global_span[0], endtime=self.global_span[1],
                timeout=eia_timeout, includerestricted=False )
        except:
            self.errmsg( "update of inventory failed, routing service failed" )
            return None
        if self.number_of_networks(slist) < eia_min_num_networks:
            self.errmsg( "update of inventory failed, number of networks %d"
                % self.number_of_networks(slist) )
            return None
        elif self._servers_missing(slist):
            self.errmsg( "update of inventory failed, servers missing %s"
                % ','.join(self._servers_missing(slist)) )
            return None
        fp = open( self.slist_cache, 'wb' )
        pickle.dump( slist, fp )
        fp.close()
        return slist
    
    def _get_inventory_from_cache( self, overrideage=False ):
        "Read station inventory from file cache."
        if not os.path.exists(self.slist_cache):
            return None
        fileage = time.time() - os.stat(self.slist_cache).st_mtime
        if fileage > self.maxcacheage and not overrideage:
            return None
        fp = open( self.slist_cache, 'rb' )
        slist = pickle.load( fp )
        fp.close()
        return slist
    
    def get_inventory( self ):
        "Read inventory from cache or from routing client."
        if self.trymgr.new_retry():
            slist = self._get_inventory_from_cache()
        else:
            slist = self._get_inventory_from_cache( overrideage=True )
        if slist is None:
            newinv = self._get_inventory_from_service()
            if newinv is None:
                self.trymgr.try_failed()
                return self._get_inventory_from_cache( overrideage=True )
            else:
                return newinv
        return slist
    
    def _servers_missing( self, inv ):
        "Check inventory for reference networks = main network of each server."
        rnets = inv.get_contents()['networks']
        miss = []
        for net in reference_networks:
            if net not in rnets:
                miss.append( net )
        return miss
    
    def number_of_networks( self, inv ):
        return len(set(inv.get_contents()['networks']))
    
    def select_random_station( self ):
        "Select random station from inventory."
        slist = self.get_inventory()
        stalist = list( set(slist.get_contents()['stations']) )
        while True:
            try:
                selsta = stalist[np.random.randint(0,len(stalist))].split()[0]
                net, sta = selsta.split('.')
            except:
                continue
            if net in large_networks.keys():
                # Throw dice to scale down probability
                if np.random.random() > large_networks[net]:
                    continue
            # Accept only operating stations in networks not excluded.
            if net not in exclude_networks and self.is_operating(slist,net,sta):
                break
        return selsta
    
    def is_operating( self, fullinv, network, station ):
        selinv = fullinv.select( network=network, station=station )
        if len(selinv) < 1:
            return False
        now = UTCDateTime()
        for episode in selinv[0]:
            if episode.end_date is None or episode.end_date > now:
                return True
        return False
    
    def _get_random_request_length( self ):
        randspan = self.maxreqlen - self.minreqlen
        return self.minreqlen + np.random.randint(0,randspan)
    
    def _get_random_request_interval( self ):
        reqspan = self._get_random_request_length()
        totalspan = int( self.global_span[1] - self.global_span[0] )
        totalspan -= reqspan
        randstart = self.global_span[0] + np.random.randint(0,totalspan)
        randend = randstart + reqspan
        return (randstart,randend)
    
    def get_station_meta( self, netsta, reqspan ):
        "Retrieve full inventory of selected station."
        net, sta = netsta.split('.')
        try:
            stamp = time.time()
            inv = self.roc.get_stations( level='response', network=net,
                station=sta, starttime=reqspan[0], endtime=reqspan[1],
                timeout=eia_timeout )
            self.meta_time = time.time() - stamp
        except Exception as e:
            self.errmsg( "requesting inventory for %s failed with %s"
                % (netsta,repr(e)))
            self.status = STATUS_NOSERV
            self.logresult( e, netsta, reqspan )
            return None
        return inv
    
    def select_random_station_channel( self, stainv, infotext='' ):
        "Select a random channel from selected station."
        wchan = [c[0:2] for c in self.wanted_channels]
        sellist = []
        for chan in set(stainv.get_contents()['channels']):
            curchan = chan.split('.')[-1]
            if curchan[0:2] in wchan:
                sellist.append( chan )
        if len(sellist) == 0:
            cont = stainv.get_contents()
            self.errmsg( "channel selection failed (%s): %s"
                % (infotext,repr(cont['channels'])) )
            return None
        return sellist[np.random.randint(0,len(sellist))]
    
    def random_request( self ):
        "Create random request parameters and return them."
        sta = self.select_random_station()
        reqspan = self._get_random_request_interval()
        stainv = self.get_station_meta( sta, reqspan )
        if stainv is None:
            return None  # Metadata request failed, error set above
        if stainv.get_contents()['channels'] == []:
            return None  # Station closed, no error
        infotext = "%s-%s" % (sta,reqspan[0].strftime('%y%m%d'))
        selchan = self.select_random_station_channel( stainv, infotext )
        if selchan is None:
            return None
        self.requestpar = (sta, selchan, stainv, reqspan)
        return (sta,selchan,stainv,reqspan)
    
    def process_request( self, station, channel, stainv, reqspan ):
        "Retrieve and evaluate waveform data, result is a status code."
        net, sta, loc, chan = channel.split('.')
        try:
            stamp = time.time()
            st = self.roc.get_waveforms( network=net, station=sta, location=loc,
                channel=chan, starttime=reqspan[0], endtime=reqspan[1] )
            self.wave_time = time.time() - stamp
        except Exception as e:
            self.errmsg( "requesting waveforms for %s failed with %s"
                % (channel,repr(e)))
            st = None
        if st is None:
            status = STATUS_NODATA
        else:
            if len(st) < 1:
                status = STATUS_NODATA
            elif len(st) > 1:
                status = STATUS_FRAGMENT
            else:
                trc = st[0]
                reqlen = reqspan[1] - reqspan[0]
                trc.trim( *reqspan )
                datalen = trc.stats.delta * len(trc.data)
                if (datalen/reqlen) < 0.99:  # FR&G&RD deliver only 0.99, why?
                    status = STATUS_INCOMPLETE
                    self.errmsg( "incomplete: %s (data:%4.2f,req:%4.2f)" % (
                        channel,datalen,reqlen) )
                else:
                    try:
                        trc.remove_response( inventory=stainv )
                        restfailed = False
                    except:
                        restfailed = True
                    if restfailed:
                        status = STATUS_RESTFAIL
                    else:
                        if (trc.data[0] != trc.data[0]):
                            status = STATUS_METAFAIL
                        else:
                            status = STATUS_OK
        self.status = status
        self.logresult()
        return (status,self.meta_time,self.wave_time)
    
    def logresult( self, exc=None, sta=None, reqspan=None ):
        "Write result of request into a file database."
        if self.requestpar is None:
            netsta = sta
            try:
                reqstart, reqend = reqspan
            except:
                reqstart = None
                reqend = None
            channel = 'unknown'
        else:
            netsta = self.requestpar[0]
            reqstart, reqend = self.requestpar[-1]
            channel = self.requestpar[1]
        net, sta = netsta.split('.')
        outpath = os.path.join( eia_datapath, 'log', net, sta )
        if not os.path.exists(outpath):
            os.makedirs( outpath )
        year = datetime.datetime.now().year
        outfile = os.path.join( outpath, "%d_%s.dat" % (year,netsta) )
        ctimestr = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        if reqstart is None:
            rtimestr = '--------_----'
            reqlen = 0.
        else:
            rtimestr = reqstart.strftime("%Y%m%d_%H%M")
            reqlen = (reqend - reqstart) / 60.
        if self.meta_time is None:
            mtimestr = "    -"
        else:
            mtimestr = "%5.1f" % self.meta_time
        if self.wave_time is None:
            wtimestr = "    -"
        else:
            wtimestr = "%5.1f" % self.wave_time
        fp = open( outfile, 'a' )
        if exc is None:
            fp.write( "%s %-8s %s %5.2f %-15s %s %s\n" % (ctimestr,
                error_names[self.status],rtimestr,reqlen,channel,mtimestr,
                wtimestr) )
        else:
            fp.write( "%s %-8s %s %5.2f %-15s %s\n" % (ctimestr,
                error_names[self.status],rtimestr,reqlen,channel,repr(exc)) )
        fp.close()
    
    def errmsg( self, text ):
        timestr = datetime.datetime.now().strftime("%d-%b-%Y_%T")
        print( "%s %s" % (timestr,text) )


#-------------------------------------------------------------------------------


class AvailabilityReport:

    # Geometry of map plot.
    mapgeometry = {
        'projection':   'merc',
        'llcrnrlon':    -11.,
        'llcrnrlat':    29.,
        'urcrnrlon':    50.,
        'urcrnrlat':    71.5,
        'resolution':   'i',
        'lat_ts':       45.,
    }
    
    # Status codes evaluated for color translation in map plot.
    okwords = ('OK',)
    failwords = ('NODATA','FRAGMENT','INCOMPL','METAFAIL','RESTFAIL')
    # Here is the file database.
    fileroot = os.path.join( eia_datapath, 'log' )

    def __init__( self ):
        self.eia = EidaAvailability()
        self.inv = self.eia.get_inventory()
        print( "inventory: "
            +"found %d networks, %d stations (with excluded networks)" % (
            len(set(self.inv.get_contents()['networks'])),
            len(set(self.inv.get_contents()['stations']))) )
        self.linecnt = 0
        self.reqstat = {}
        self.netstat = {}
        self.current_network = None
        self.repfp = None
        self.etime = datetime.datetime.now()
        self.etimestr = self.etime.strftime("%d-%b-%Y_%H:%M")
        self.stime = self.etime \
            - datetime.timedelta( days=eia_reqstats_timespan_days )
        self.minreqtime = None
        self.report_outpath = None
    
    def parse_yearfile( self, fname, okcnt, failcnt ):
        """Parse a single file in the file database and return the sum of
        evaluated status codes as well as a location.
        """
        lat = lon = None
        for line in open(fname):
            tmp = line.split()
            if len(tmp) < 2:
                continue
            try:
                reqtime = datetime.datetime.strptime( tmp[0], "%Y%m%d_%H%M" )
            except:
                print( "Error parsing file '%s'" % fname )
                continue
            if reqtime < self.stime:
                # If waveform request is too old, ignore.
                continue
            if self.minreqtime is None or self.minreqtime > reqtime:
                self.minreqtime = reqtime
            keyw = tmp[1]
            if lat is None and len(tmp) > 4:
                chan = tmp[4]
                chan = chan[:-1] + 'Z'
                try:
                    coo = self.inv.get_coordinates( chan )
                except:
                    #if not chan.startswith('unknow'):
                    #    pass
                    continue
                if 'latitude' in coo.keys():
                    lat = coo['latitude']
                if 'longitude' in coo.keys():
                    lon = coo['longitude']
            if keyw in self.okwords:
                okcnt += 1
            elif keyw in self.failwords:
                failcnt += 1
            self.linecnt += 1
            self.add_keyword( keyw )
        return (okcnt,failcnt,lat,lon)

    def add_keyword( self, keyw ):
        "Store status codes for network statistics."
        if keyw not in error_names.values():
            return
        if not self.current_network in self.netstat.keys():
            self.netstat[self.current_network] = {}
            for errmsg in error_names.values():
                self.netstat[self.current_network][errmsg] = 0
        self.netstat[self.current_network][keyw] += 1

    def parse_years( self, fpath ):
        "Parse all files of a station."
        okcnt = failcnt = 0
        lincnt = self.linecnt
        for fname in sorted(os.listdir(fpath)):
            fullname = os.path.join( fpath, fname )
            okcnt, failcnt, lat, lon = self.parse_yearfile( fullname,
                okcnt, failcnt )
        if okcnt == 0 and failcnt == 0:
            return (None,None,None)
        okperc = 100. * float(okcnt) / float(okcnt+failcnt)
        self.add_stats( self.linecnt - lincnt )
        return (okperc,lat,lon)
    
    def loop_files( self ):
        """Loop all networks and stations in file database, return availability
        and location.
        """
        data = []
        for netdir in sorted(os.listdir(self.fileroot)):
            if len(netdir) != 2:
                continue
            self.current_network = netdir
            stapath = os.path.join( self.fileroot, netdir )
            for stadir in sorted(os.listdir(stapath)):
                fpath = os.path.join( stapath, stadir )
                okperc, lat, lon = self.parse_years( fpath )
                if None in (okperc,lat,lon):
                    continue
                data.append( (okperc,lat,lon) )
        return data
    
    def add_stats( self, cnt ):
        "Hit count statistics, how many stations have how many hits."
        if not cnt in self.reqstat.keys():
            self.reqstat[cnt] = 0
        self.reqstat[cnt] += 1
    
    def total_number_of_stations( self ):
        """Determine total number of stations, remove double entries and
        stations of excluded networks.
        """
        stalist = []
        double_entries = 0
        for statext in self.inv.get_contents()['stations']:
            netsta = statext.split()[0]
            net, sta = netsta.split('.')
            if net in exclude_networks:
                continue
            if netsta in stalist:
                double_entries += 1
            else:
                stalist.append( netsta )
        return (len(stalist),double_entries)
    
    def dump_netstat( self ):
        "Write networks status statistics in different formats."
        def print_header( skeys ):
            xkeys = ["`%s`" % k for k in skeys]
            header = '|net ' + ' '.join(["| %10s" % k for k in xkeys]) + ' |'
            self.repprint( '+----+' + '+'.join(len(skeys)*['------------']) + '+' )
            self.repprint( header )
            self.repprint( '+:===+' + '+'.join(len(skeys)*['===========:']) + '+' )
        skeys = sorted( error_names.values() )
        skeys.remove( error_names[STATUS_OK] )
        skeys.remove( error_names[STATUS_NODATA] )
        skeys = [error_names[STATUS_OK],error_names[STATUS_NODATA]] + skeys
        self.newpage()
        tableheader = "Request status statistics of networks"
        maxlen = 45
        for netnum,net in enumerate(sorted(self.netstat.keys())):
            if netnum % maxlen == 0:
                self.repprint( "\n" )
                self.repprint( "\n**%s:**\n" % tableheader )
                if netnum == 0:
                    tableheader += " (continued)"
                print_header( skeys )
            netline = "| %s " % net \
                + ' '.join(["| %10d" % self.netstat[net][k] for k in skeys]) + ' |'
            self.repprint( netline )
        self.repprint( '+----+' + '+'.join(len(skeys)*['----------']) + '+' )
        self.repprint( "\nStatus codes used in above statistics:\n" )
        self.repprint( "`OK`       \n: data delivery and restitution successful\n" )
        self.repprint( "`NODATA`   \n: no data available\n" )
        self.repprint( "`FRAGMENT` \n: returned data not contigous\n" )
        self.repprint( "`INCOMPL`  \n: returned time interval less than requested\n" )
        self.repprint( "`METAFAIL` \n: restituted data contain illegal values (`Nan`s)\n" )
        self.repprint( "`NOSERV`   \n: station metadata request failed\n" )
        self.repprint( "`RESTFAIL` \n: removing response failed\n" )
    
    def makeavailplot( self, outfile=None, mapgeo=None ):
        plt.figure( figsize=(14,10) )
        if mapgeo is None:
            xmap = Basemap( **self.mapgeometry )
        else:
            xmap = Basemap( **mapgeo )
        xmap.drawcoastlines(linewidth=0.25, zorder=3)
        xmap.drawcountries(linewidth=0.25, zorder=3)
        xmap.fillcontinents( color='#FFFFFF', lake_color='#EEEEFF', zorder=2 )
        xmap.drawmapboundary( fill_color='#EEFEFF' )
        data = self.loop_files()
        ##cols, lats, lons = zip( *data )  # plot nice colors at last.
        # Plot first the red dots, then overlain by other colors.
        redcols, redlats, redlons = ([],[],[])
        nonredcols, nonredlats, nonredlons = ([],[],[])
        for col,lat,lon in data:
            if col < 0.1:
                redcols.append( col )
                redlats.append( lat )
                redlons.append( lon )
            else:
                nonredcols.append( col )
                nonredlats.append( lat )
                nonredlons.append( lon )
        xv, yv = xmap( redlons, redlats )
        xmap.scatter( xv, yv, c=redcols, edgecolor='', cmap='RdYlGn', s=10,
            zorder=5 )
        xv, yv = xmap( nonredlons, nonredlats )
        xmap.scatter( xv, yv, c=nonredcols, edgecolor='', cmap='RdYlGn', s=10,
            zorder=6 )
        datestr = datetime.datetime.now().strftime("%y%m%d")
        plt.title( "EIDA waveform response statistics (%s)" % datestr )
        plt.draw()
        if outfile:
            plt.savefig( outfile )
            tmpfile = os.path.join( os.path.dirname(outfile),
                'xxx_'+os.path.basename(outfile) )
            shellcmd = "convert %s -trim %s; mv %s %s" % (outfile,
                tmpfile,tmpfile,outfile)
            #print( "dbg: shellcmd", shellcmd )
            os.system( shellcmd )
        else:
            plt.show()
        return len(data)
    
    def makehitplot( self, hitstats, outfile=None ):
        fig = plt.figure( figsize=(8,6) )
        ax = fig.add_subplot( 111 )
        xvals = sorted( hitstats.keys() )
        yvals = [hitstats[x] for x in xvals]
        ax.plot( xvals, yvals, c='r' )
        ax.set_xlabel( "number of hits" )
        ax.set_ylabel( "number of stations" )
        if outfile:
            plt.savefig( outfile )
            tmpfile = os.path.join( os.path.dirname(outfile),
                'xxx_'+os.path.basename(outfile) )
            shellcmd = "convert %s -trim %s; mv %s %s" % (outfile,
                tmpfile,tmpfile,outfile)
            #print( "dbg: shellcmd", shellcmd )
            os.system( shellcmd )
        else:
            plt.show()
    
    def make_md_report( self, repfile ):
        "Create pandoc markdown report file."
        def add_inventory_test_report( fp, sday, respfig ):
            """Retrieve statistics on the EIDA inventory test and include it
            in the report. The code 'eida_inventory_test.py' is run on a
            different host and we need access to its log file. Therefore this
            part contains implementation specific code.
            """
            remhost = eia_spec_invtest_host
            rempath = eia_spec_invtest_rempath
            logfilename = "eida_inventory_test.log"
            shellcmd = "scp -q %s:%s/%s ." % (remhost,rempath,logfilename)
            ret = os.system( shellcmd )
            if ret != 0:
                # No access to logfile of eida_inventoy_test.py, omit this part.
                print( "*** eida_inventory_test logfile not found. ***" )
                return None
            # The parser program is on the same path as this one.
            parseprog = os.path.join( os.path.dirname(sys.argv[0]),
                'parse_eida_inventory_test_log.py' )
            shellcmd = "%s %s %s report %s" % (parseprog,logfilename,
                sday.strftime("%d-%b-%Y_%T"),respfig)
            #print( "dbg: executing:", shellcmd )
            invtest_data = ""
            for line in os.popen( shellcmd ).readlines():
                invtest_data += line
            if len(invtest_data) < 10:
                # Should not happen, something with the parsing went wrong.
                print( "*** Could not parse inventory test logfile. ***" )
                if os.path.exists(logfilename):
                    os.remove( logfilename )
                return None
            invtest_text = "## Failure rate of inventory requests\n\n" \
                + "This section contains results of inventory test \n" \
                + "requests on network, station and channel level.\n" \
                + "A few times per hour all servers get direct\n" \
                + "metadata requests followed by a metadata request\n" \
                + "using the routing client of obspy. It is checked\n" \
                + "whether all servers respond to the direct requests\n" \
                + "and whether all servers contribute to the routed\n" \
                + "request. The following results refer to tests carried\n" \
                + "out since %s.\n\n" % sday.strftime("%d-%b-%Y %T")
            respfigtext = "\n![Responsitivity of all servers plotted with a " \
                +"granularity of 8h; green = 0% errors, orange = 10%, " \
                +"brown = 50%%, black = 100%%](%s){width=100%%}" % (
                os.path.basename(respfig))
            if os.path.exists(logfilename):
                os.remove( logfilename )
            return invtest_text+invtest_data+respfigtext
        figfile = os.path.splitext(repfile)[0] + '_figure.png'
        numstations = self.makeavailplot( figfile )
        figfile2 = os.path.splitext(repfile)[0] + '_fig2.png'
        figfile3 = os.path.splitext(repfile)[0] + '_fig3.png'
        self.makehitplot( self.reqstat, figfile2 )
        begtime = max( self.minreqtime, self.stime )
        begtimestr = begtime.strftime( "%d-%b-%Y" )
        fp = open( repfile,  'w' )
        self.repfp = fp
        self.repprint( "# EIDA Availability Report\n" )
        self.repprint( "## Created at %s\n" % self.etimestr.replace('_',' ') )
        self.repprint( "This document contains results of automated tests\n"
            +"of the waveform availability of European EIDA stations and the\n"
            +"responsitivity of the EIDA servers to metadata requests.\n"
            )

        self.repprint( "## Description of waveform test program" )
        self.repprint( __doc__ )
        self.repprint( "\n## Statistics on waveform tests\n" )
        self.repprint( "Statistics on random requests between %s and %s" % (
            begtimestr,self.etimestr.replace('_',' ')) )
        self.repprint( "using station metadata valid since %s.\n" % 
            (datetime.datetime.now()-datetime.timedelta(days=365)).strftime("%d-%b-%Y") )
        totstations, doublesta = self.total_number_of_stations()
        self.repprint( "\n**Counters:**\n" )
        self.repprint( "- unrestricted stations offering channels `%s`: %d"
            % (','.join(self.eia.wanted_channels),totstations) )
        self.repprint( "- evaluated stations: %d" % numstations )
        self.repprint( "- number of requests: %d" % self.linecnt )
        self.newpage()

        self.repprint( "## Waveform availability plot\n" )
        self.repprint( "Color coded plot of evaluated EIDA stations. " )
        self.repprint( "Shows results of %d random requests between %s and %s."
            % (self.linecnt,begtimestr,self.etimestr.split('_')[0]) )
        self.repprint( "The availability displayed is computed as the relative number\n"
            +"of request results with status `OK` (see table below) compared\n"
            +"to the number of all requests to this station.\n"
        )
        self.repprint( "![Availability of stations: green 100%, yellow 50%, "\
            +"red 0%%](%s){width=100%%}" % os.path.basename(figfile) )
        self.newpage()
        
        #self.repprint( "## Network statistics\n" )
        self.dump_netstat()
        self.newpage()

        self.repprint( "\n## Waveform requests, random hit distribution\n" )
        self.repprint( "How many stations have how many hits of random " )
        self.repprint( "requests.\n" )
        #self.repprint( "  No. of hits   Stations with this number of hits" )
        #self.repprint( "-------------   -----------------------------------------------" )
        #for k in sorted(self.reqstat.keys()):
        #    self.repprint( "       %4d    %4d" % (k,self.reqstat[k]) )
        self.repprint( "![Request hit statistics showing the distribution of "
            +"the %d requests on the %d evaluated stations](%s){width=100%%}"
            % (self.linecnt,numstations,os.path.basename(figfile2)) )
        self.repprint( "\n" )
        self.newpage()

        sday = datetime.datetime.now() - datetime.timedelta( days=30 )
        invtest_report = add_inventory_test_report( fp, sday, figfile3 )
        if invtest_report:
            self.newpage()
            self.repprint( invtest_report )

        self.repprint( "\n\n## Remarks\n\n" )
        self.repprint( "A history of these daily reports (in pdf format)" )
        self.repprint( "as well as request logs on station level are available at " )
        self.repprint( "<ftp://www.szgrf.bgr.de/pub/EidaAvailability>," )
        self.repprint( "files `history_eida_availability_reports.tgz` and " )
        self.repprint( "`stationlogs_eida_availability.tgz`, respectively." )
        self.repprint( "\n\nThis report was automatically created at %s MEST using"
            % self.etimestr.replace('_',' ') )
        self.repprint( "%s.\n" % os.popen('pandoc --version').readline().strip() )
        fp.close()
        self.reffp = None
    
    def make_pdf_report( self, mdreport ):
        pdffile = os.path.splitext(mdreport)[0] + '.pdf'
        shellcmd = "cd %s; pandoc -V papersize=a4 " % os.path.dirname(mdreport)\
            +"-V mainfont=Verdana "\
            +"-V fontsize=10pt --pdf-engine xelatex "\
            +"-o %s %s" % (pdffile,mdreport)
        # alternative fonts: Verdana, NimbusSanL-Regu
        #print ("dbg: shellcmd", shellcmd )
        print( "Creating pdf file '%s'" % pdffile )
        os.system( shellcmd )
        return pdffile
    
    def make_html_report( self, mdreport, cssfile=None ):
        if cssfile is None:
            cssfile = eia_spec_default_cssfile
        if not os.path.exists(cssfile):
            print( "Need existing css file for HTML output" )
            print( "Specified file '%s' not found" % cssfile )
            return
        if self.report_outpath:
            os.system( "cp %s %s/" % (cssfile,self.report_outpath) )
        else:
            os.system( "cp %s ." % cssfile )
        lcssfile = os.path.basename( cssfile )
        htmlfile = os.path.splitext(mdreport)[0] + '.html'
        htmltitle = "EIDA Availability Report"
        shellcmd = "cd %s; pandoc -s -c %s --metadata title='%s' -o %s %s" % (
            os.path.dirname(mdreport),lcssfile,htmltitle,htmlfile,mdreport)
        #print ("dbg: shellcmd", shellcmd )
        print( "Creating HTML file '%s'" % htmlfile )
        os.system( shellcmd )
        return htmlfile
    
    def display_pdf_report( self, pdffile ):
        os.system( "evince %s" % pdffile )
    
    def repprint( self, text ):
        if self.repfp is None:
            print( text )
        else:
            self.repfp.write( "%s\n" % text )
    
    def newpage( self ):
        self.repprint( "\n\n\\newpage\n\n" )
    
    def daily_report( self, outpath=None ):
        self.report_outpath = outpath
        repname = "eida_availability_report.md"
        if outpath:
            frepname = os.path.join( outpath, repname )
        else:
            frepname = repname
        self.make_md_report( frepname )
        self.make_pdf_report( frepname )
        self.make_html_report( frepname )
        #os.remove( frepname )
        


#-------------------------------------------------------------------------------


class DoubleProcessCheck:

    def __init__( self ):
        self.pidfile = '/tmp/EidaAvailability.pid'
        self.maxage = 300
    
    def create_pidfile( self ):
        fp = open( self.pidfile, 'w' )
        fp.write( "%d\n" % os.getpid() )
        fp.close()
    
    def process_active( self ):
        if not os.path.exists(self.pidfile):
            return False
        fileage = time.time() - os.stat(self.pidfile).st_mtime
        if fileage < self.maxage:
            return True
        try:
            pid = int( open(self.pidfile).readline().strip() )
        except:
            os.remove( self.pidfile )
            return False
        try:
            timestr = datetime.datetime.now().strftime("%d-%b-%Y_%T")
            print( "DoubleProcessCheck: %s killing process %d" % (timestr,pid) )
            os.kill( pid, signal.SIGTERM )
            time.sleep( 5 )
            os.kill( pid, signal.SIGKILL )
        except:
            pass
        os.remove( self.pidfile )
        return False
    
    def should_exit( self ):
        if self.process_active():
            return True
        self.create_pidfile()
        return False
    
    def release( self ):
        if os.path.exists(self.pidfile):
            os.remove(self.pidfile)


#-------------------------------------------------------------------------------


class RetryManager:

    def __init__( self, name, waittime=3600 ):
        self.name = name
        self.waittime = float(waittime)
        self.flagfile = "/tmp/retrymanager_%s.flag" % self.name
    
    def try_failed( self ):
        "Mark a failed try by touching the flag file."
        os.system( "touch %s" % self.flagfile )
    
    def new_retry( self ):
        "Return True/False depending on age of flag file."
        if not os.path.exists(self.flagfile):
            return True
        fileage = time.time() - os.stat(self.flagfile).st_mtime
        retry = (fileage > self.waittime)
        if retry:
            os.remove( self.flagfile )
        return retry


#-------------------------------------------------------------------------------


if __name__ == '__main__':

    pcheck = DoubleProcessCheck()
    if pcheck.should_exit():
        exit()

    stamp = time.time()
    eia = EidaAvailability()
    rr = eia.random_request()
    if rr is not None:
        eiaresult = eia.process_request( *rr )
        runtime = time.time() - stamp
        print( "status %s %s %3.1f" % (rr[1], repr(eiaresult), runtime ) )

    pcheck.release()
