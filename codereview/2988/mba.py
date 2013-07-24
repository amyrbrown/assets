# Copyright (C) 2012 Sergio Rossell
#
# This script is part of the EXAMO software
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
# 

"""
MBA implementation
"""
from numpy import array
from gurobipy import *
import random 
import time 
import os

from sys import path; path.append('./modules/')
from utilities import importPickle, exportPickle
from examo import CbModel, deleteCbmRxns



################################################################################
# FUNCTIONS

def findActiveRxns(cbm, thresh, rl = []):
    act = set()
    arrayIdRs = array(cbm.idRs[:])
    cbm.initLp()
    if rl:
        idRs = rl
    else:
        idRs = cbm.idRs[:]
    # maximizing all reactions at once
    # reseting the objective
    cbm.guro.setObjective(0)
    # setting the objective
    s = 'cbm.linobj = LinExpr([1.0] * len(cbm.idRs), ['
    for var in cbm.guro.getVars():
        s += 'cbm.%s, ' % var.varName
    s = s.rstrip(', ')
    s += '])'
    exec s
    cbm.guro.setObjective(cbm.linobj, -1)#1 for maximize
    cbm.guro.optimize()
    sol = abs(array([v.x for v in cbm.guro.getVars()]))
    indices = (sol > thresh).nonzero()[0]
    act.update(arrayIdRs[indices])
    idRs = list(set(idRs) - act)
    # maximizing
    for rxn in idRs:
        if rxn not in act:
            # reseting the objective
            cbm.guro.setObjective(0)
            exec 'cbm.guro.setObjective(cbm.%s, GRB.MAXIMIZE)' % rxn
            cbm.guro.optimize()
            sol = abs(array([v.x for v in cbm.guro.getVars()]))
            indices = (sol > thresh).nonzero()[0]
            act.update(arrayIdRs[indices])
    idRs = list(set(idRs) - act)
    # minimizing
    for rxn in idRs:
        if rxn not in act:
            # reseting the objective
            cbm.guro.setObjective(0)
            exec 'cbm.guro.setObjective(cbm.%s, GRB.MINIMIZE)' % rxn
            cbm.guro.optimize()
            sol = abs(array([v.x for v in cbm.guro.getVars()]))
            indices = (sol > thresh).nonzero()[0]
            act.update(arrayIdRs[indices])
    return act

def pruneRxn(cbm, cH, rxn, thresh):
    m0 = deleteCbmRxns(cbm, rxn)
    #NOTE the threshold for is set a bit higher for cH rxns
    act = findActiveRxns(m0, thresh, cH)
    if len(cH - act) != 0:#not all cH rxns are active
        #print len(cH - act)
        return cbm
    else:
        act.update(findActiveRxns(m0, thresh))
        inact = set(m0.idRs) - act
        inact.add(rxn)
        m1 = deleteCbmRxns(m0, inact)
        return m1


def iterativePrunning(m, cH, thresh = 1E-10, prunableRxns = set()):
    """
    solver can be 'cplex', 'glpk' or 'gurobi'
    """
    if prunableRxns:
        prunableRxns = list(set(prunableRxns) - cH)
    else:
        prunableRxns = list(set(m.idRs) - cH)
    semilla = int((time.time() * 1E6) * os.getpid())
    random.seed(semilla)
    random.shuffle(prunableRxns)
    prunableRxns = set(prunableRxns)
    while prunableRxns:
        rxn = prunableRxns.pop()
        try:
            mTemp = pruneRxn(mTemp, cH, rxn, thresh)
            prunableRxns = prunableRxns & set(mTemp.idRs)
        except NameError:
            mTemp = pruneRxn(m, cH, rxn, thresh)
            prunableRxns = prunableRxns & set(mTemp.idRs)
        #print os.getpid(), len(prunableRxns)
    return mTemp.idRs

if __name__ == '__main__':
    ########################################
    # INPUTS
    
    md = importPickle('../data/iMM904_blkRxnsDeleted_dict.pkl')
    rH = importPickle('../data/rxnsClassifiedByExprssion_eth_15_85.pkl')['rH']
    fOutMbaCandRxns = '../data/mbaCandRxns/mbaCandRxns_%s.pkl'

    activityThreshold = 1E-10
    descrip = 'eth_15_85'
    ########################################
    # STATEMENTS
    m = CbModel(md['S'], md['idSp'], md['idRs'], md['lb'], md['ub'], md['rxns'],
        md['genes'])
    #act = findActiveRxns(m, activityThreshold, m.idRs)
    #mT = pruneRxn(m, rH, 'R_PFK', activityThreshold)
    #iterativePrunning(m, rH, activityThreshold)
    
    numProc = 2
    numRep = 2
    # making sure that all rH reactions are active to begin with
    act = findActiveRxns(m, 1E-10, rH)
    cH = rH & act
    cH.add('R_biomass_published')
    print '%i processes and %i repetitions' % (numProc, numRep)
    for i in range(numProc):
        pid = os.fork()
        if pid == 0:
            for j in range(numRep):
                locTime = time.localtime()
                pid = os.getpid()
                timeStr = '%i%02i%02i%02i%02i%02i' % locTime[:6]
                tag = '%s_%s_%s' % (descrip, pid, timeStr)
                cr = iterativePrunning(m, cH, activityThreshold, m.idRs[:10])
                exportPickle(cr, fOutMbaCandRxns % tag)
            os._exit(0)




