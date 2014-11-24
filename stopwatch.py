#!/usr/bin/env python
"""
StopWatch
========

Created on Thu Sep 27 11:38:10 2012
by fwagner@gfz-potsdam.de
"""

import time


def start(string="Starting computation..."):
    print(string)
    global t_0
    t_0 = time.time()


def stop():
    elapsed = time.time() - t_0
    m, s = divmod(elapsed, 60)
    h, m = divmod(m, 60)
    print(("Elapsed run time: %02dh %02dmin %02dsec" % (h, m, s)))
    return elapsed
