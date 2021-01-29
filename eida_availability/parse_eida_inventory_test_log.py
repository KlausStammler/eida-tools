#! /usr/bin/env python
#
# file parse_eida_inventory_test_log.py
#      ================================
#
# K. Stammler, 26-Jun-2020

"""
Parse logfile of eida_inventory_test.py

Syntax:
    parse_eida_inventory_test_log.py <logfile> [<starttime>] [<mode>] [<plotname>]
    
    <logfile>   logfile of program 'eida_inventory_test.py'
    <starttime> start time of evaluation, default: beginning of file
                format example: 23-Jun-2020_15:00
    <mode>      'normal' (default) or 'report'
    <plotname>  if not empty create plot with response info per server

Example:
    parse_eida_inventory_test_log.py ~/Work/EidaAvailability/eida_inventory_test.log
"""

from __future__ import print_function
import os
import sys
import datetime
if not 'DISPLAY' in os.environ.keys():
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

reference_networks = {
    'NL' : 'ODC',
    'GE' : 'GFZ',
    'FR' : 'RESIF',
    'CH' : 'ETH',
    'GR' : 'BGR',
    'BW' : 'LMU',
    'RO' : 'NIEP',
    'KO' : 'KOERI',
    'HL' : 'NOA',
    'NO' : 'UIB/NORSAR',
    'CA' : 'ICGC',
    'IV' : 'INGV',
}


#-------------------------------------------------------------------------------


def parse_logfile( logfile, stime=None, mode='normal', respplot=None ):

    skip = False
    roclifailures = 0
    noerrorcnt = 0
    directcnt = 0
    routecnt = 1  # first is missed in loop
    directfail = {}
    routefail = {}
    routeactive = False
    currtime = None
    currsrv = None
    missnet = ""
    failedlist = []
    if stime and respplot:
        mrp = MetaResponsePlot( stime, 8, respplot )
    else:
        mrp = None
    for line in open(logfile):
        if line.startswith('eida_inventory_test.py started at'):
            timestr = line.split()[3]
            currtime = datetime.datetime.strptime( timestr, '%d-%b-%Y_%H:%M:%S' )
            if stime is not None and currtime < stime:
                skip = True
                continue
            else:
                skip = False
            if mrp:
                mrp.set_time( currtime )
            directcnt += 1
            if routeactive:
                routecnt += 1
                if mrp:
                    mrp.inc_total_count( 'ALL' )
        elif skip:
            continue
        elif 'reading inventory from server' in line:
            srv = line.split()[4]
            if srv == 'http://eida.geo.uib.no':
                srv = 'UIB/NORSAR'
            elif srv == 'https://eida.bgr.de':
                srv = 'BGR'
            if mrp:
                mrp.inc_total_count( srv )
        elif 'reading inventory from routing client' in line:
            srv = None
        elif 'FAILED:' in line:
            if srv is None:
                roclifailures += 1
            else:
                failedlist.append( srv )
                if mrp:
                    mrp.inc_fail_count( srv )
        elif line.startswith('missing reference networks'):
            routeactive = True
            missnet = line.split()[3]
        elif line.startswith('============'):
            timestr = currtime.strftime( "%d-%b-%Y_%T" )
            tmissnet = transref(missnet)
            if not failedlist and not missnet:
                noerrorcnt += 1
            if mode == 'normal':
                print( "%s %20s %20s" % (timestr,','.join(failedlist),tmissnet) )
            cumulate_failures( failedlist, tmissnet.split(','), directfail,
                routefail )
            if tmissnet and mrp:
                for srv in tmissnet.split(','):
                    mrp.inc_fail_count( srv )
            failedlist = []
            missnet = ""
    
    if mrp: mrp.finish_index()
    print( "totals: direct requests %d, routed requests %d" % (directcnt,routecnt) )
    print_failure_rates( directfail, routefail, directcnt, routecnt,
        noerrorcnt, roclifailures, mode )
    if mrp:
        mrp.makeplot()
        #mrp.dump()


#-------------------------------------------------------------------------------


def transref( missnet ):
    if missnet == '':
        return ''
    srvlist = []
    for net in missnet.split(','):
        srvlist.append( reference_networks[net] )
    #if 'ICGC' in srvlist:
    #    srvlist.remove( 'ICGC' )
    return ','.join(sorted(srvlist))


#-------------------------------------------------------------------------------


def print_failure_rates( directfail, routefail, directcnt, routecnt,
    noerrcnt, roclifailures, mode ):
    print( "\nNumber of failed requests and failure rates of servers:" )
    print()
    if mode == 'normal':
        print( "server    direct          routed" )
    else:
        print()
        print( "| server  |  direct       |  routed       |" )
        print( "|--------:|--------------:|--------------:|" )
    for srv in sorted(reference_networks.values()):
        #if srv == 'ICGC': continue
        dnum = directfail.get( srv, 0 )
        dperc = 100.0*float(dnum)/float(directcnt)
        rnum = routefail.get( srv, 0 )
        rperc = 100.0*float(rnum)/float(routecnt)
        if mode == 'normal':
            line = " %5s %5d (%4.1f%%)  %5d (%4.1f%%)" % (srv,dnum,dperc,rnum,rperc)
            print( line )
        else:
            line = "| %5s   | %5d (%4.1f%%) | %5d (%4.1f%%) |" % (srv,dnum,dperc,rnum,rperc)
            print( line )
    noerrperc = 100.*float(noerrcnt)/float(directcnt)
    print( "\nfailures of routing client: %d" % roclifailures )
    if mode == 'report':
        print()
    print( "runs without errors: %d (%3.1f%%)" % (noerrcnt,noerrperc) )

#-------------------------------------------------------------------------------


def cumulate_failures( dfail, rfail, ddict, rdict ):
    for srv in dfail:
        if not srv in ddict.keys():
            ddict[srv] = 0
        ddict[srv] += 1
    for srv in rfail:
        if not srv in rdict.keys():
            rdict[srv] = 0
        rdict[srv] += 1


#-------------------------------------------------------------------------------


def parse_time( timestr ):
    try:
        stime = datetime.datetime.strptime( timestr, "%d-%b-%Y_%H:%M" )
    except:
        stime = None
    if stime is None:
        try:
            stime = datetime.datetime.strptime( timestr, "%d-%b-%Y_%H:%M:%S" )
        except:
            stime = None
    if stime is None:
        try:
            stime = datetime.datetime.strptime( timestr, "%d-%b-%Y" )
        except:
            stime = None
    return stime


#-------------------------------------------------------------------------------


class MetaResponsePlot:

    """Create a plot with response info of all servers."""

    def __init__( self, stime, granularity, plotname ):
        self.stime = stime   # start time of plot, end time is now
        self.granularity = float(granularity)   # in hours
        self.plotname = plotname
        self.totcnt = {}   # counter per server
        self.failcnt = {}  # counter per server
        self.currtime = None
        self.curridx = None
        self.lastidx = None
        self.plotinfo = {}  # dict of dicts
        self.incall_on_none = False
    
    def set_time( self, ctime ):
        if ctime == self.currtime:
            return
        self.currtime = ctime
        difftime = ctime - self.stime
        diffsec = difftime.seconds + difftime.days*86400
        self.curridx = int((diffsec/3600.)/self.granularity)
        if self.lastidx != self.curridx:
            self.finish_index()
            self.lastidx = self.curridx
    
    def time_label( self, tm ):
        difftime = tm - self.stime
        diffsec = difftime.seconds + difftime.days*86400
        return ((diffsec/3600.)/self.granularity)
    
    def inc_total_count( self, server ):
        if server == 'ALL':
            if self.totcnt == {}:
                self.incall_on_none = True
            else:
                for k in self.totcnt.keys():
                    self.totcnt[k] += 1
                if self.incall_on_none:
                    for k in self.totcnt.keys():
                        self.totcnt[k] += 1
                    self.incall_on_none = False
            return
        if not server in self.totcnt.keys():
            self.totcnt[server] = 0
        self.totcnt[server] += 1

    def inc_fail_count( self, server ):
        if not server in self.failcnt.keys():
            self.failcnt[server] = 0
        self.failcnt[server] += 1
    
    def finish_index( self ):
        for k in self.totcnt.keys():
            if k in self.failcnt.keys():
                fcnt = self.failcnt[k]
            else:
                fcnt = 0
            if k not in self.plotinfo.keys():
                self.plotinfo[k] = {}
            #if k == 'NIEP':
            #    print( "dbg: fcnt, totcnt", fcnt, self.totcnt[k] )
            self.plotinfo[k][self.curridx] = float(fcnt)/float(self.totcnt[k])
        self.totcnt = {}
        self.failcnt = {}
    
    def makeplot( self ):
        fig = plt.figure( figsize=(8,6) )
        ax = fig.add_subplot( 111 )
        ypos = 0.5
        ylabs = []
        yticks = []
        for srv in sorted(self.plotinfo.keys(),reverse=True):
            for idx in range(min(self.plotinfo[srv]),max(self.plotinfo[srv])+1):
                ax.plot( (idx,idx+0.2), (ypos,ypos), c=self.getcol(srv,idx),
                    linewidth=9 )
            yticks.append( ypos )
            ylabs.append( srv )
            ypos += 0.5
        plt.yticks( yticks, ylabs )
        xtickpos = datetime.datetime( self.stime.year, self.stime.month,
            self.stime.day, 0 ) + datetime.timedelta( days=1 )
        xticks = []
        xlabs = []
        while xtickpos < datetime.datetime.now():
            xticks.append( self.time_label(xtickpos)+1. )
            xlabs.append( xtickpos.strftime("%d-%b") )
            xtickpos += datetime.timedelta( days=4 )
        ax.set_title( "responsitivity to metadata requests (%d)"
            % self.stime.year )
        plt.xticks( xticks, xlabs )
        #plt.show()
        plt.savefig( self.plotname )
    
    def getcol( self, srv, idx ):
        if not srv in self.plotinfo.keys():
            return 'gray'
        elif not idx in self.plotinfo[srv].keys():
            return 'gray'
        val = self.plotinfo[srv][idx]
        col1r = 0.
        col1g = 1.0
        col2r = 1.0
        col2g = 0.5
        col3r = 0.
        col3g = 0.
        thresh1 = 0.1
        if val < thresh1:
            fac = val/thresh1
            colr = col1r + fac * (col2r-col1r)
            colg = col1g + fac * (col2g-col1g)
        else:
            fac = (val-thresh1)/(1.-thresh1)
            colr = col2r + fac * (col3r-col2r)
            colg = col2g + fac * (col3g-col2g)
        xcolr = int(colr * 255)
        xcolg = int(colg * 255)
        return '#%02x%02x00' % (xcolr,xcolg)
    
    def dump( self ):
        print( "dbg: mrp plotinfo:", self.plotinfo )


#-------------------------------------------------------------------------------


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print( __doc__ )
        exit()
    
    argc = 1
    logfile = sys.argv[argc]
    argc += 1
    stime = None
    if len(sys.argv) > argc and sys.argv[argc]:
        stime = parse_time( sys.argv[argc] )
        if stime is None:
            print( "Cannot parse time string '%s'" % sys.argv[argc] )
            exit()
    argc += 1
    mode = 'normal'
    if len(sys.argv) > argc and sys.argv[argc]:
        mode = sys.argv[argc]
    if mode not in ('normal','report'):
        print( "Illegal mode '%s'" % mode )
        exit()
    argc += 1
    respplot = None
    if len(sys.argv) > argc and sys.argv[argc]:
        respplot = sys.argv[argc]
    
    
    if not os.path.exists(logfile):
        print( "Cannot find input file '%s'" % logfile )
        exit()
    
    parse_logfile( logfile, stime, mode, respplot )
