import numpy as np
import os, glob, json
import optparse

import enterprise
from enterprise.pulsar import Pulsar

from enterprise.signals.parameter import function
from enterprise.signals.deterministic_signals import Deterministic
from enterprise.signals import parameter, selections
from enterprise.signals.selections import Selection

from enterprise_extensions import models
from enterprise_extensions.frequentist import optimal_statistic as op

parser = optparse.OptionParser()
parser.add_option('--datadir', action='store', dest='datadir', default='./', type='string')
parser.add_option('--outdir', action='store', dest='outdir', default='./report/gwb_run/', type='string')
parser.add_option('--noisedir', action='store', dest='noisedir', default='./noisefiles/', type='string')
parser.add_option('--orf', action='store', dest='orf', default='crn', type='string')
parser.add_option('--common_psd', action='store', dest='common_psd', default='powerlaw', type='string')
parser.add_option('--common_components', action='store', dest='common_components', default=30, type='int')
parser.add_option('--gamma_common', action='store', dest='gamma_common', default=4.33)
parser.add_option('--red_components', action='store', dest='red_components', default=0, type='int')
parser.add_option('--dm_components', action='store', dest='dm_components', default=0, type='int')
parser.add_option('--chrom_components', action='store', dest='chrom_components', default=0, type='int')
parser.add_option('--num_dmdips', action='store', dest='num_dmdips', default=2, type='int')
parser.add_option('--bayesephem', action='store_true', dest='bayesephem', default=False)
parser.add_option('--common_sin', action='store_true', dest='sin_wave', default=False)
parser.add_option('--number', action='store', dest='num', default=1e4, type='float')
(options,args) = parser.parse_args()

# load par and tim files
datadir = options.datadir
outdir = options.outdir
noisedir = options.noisedir

parfiles = sorted(glob.glob(datadir + '/J*/*.par'))
timfiles = sorted(glob.glob(datadir + '/J*/*_all.tim'))
noisefiles = sorted(glob.glob(noisedir + '/*.json'))

# filter to one set of par+tim+noisefile per pulsar
PsrList = np.loadtxt(datadir + 'scripts/psrlist.txt',dtype=str)

parfiles = [x for x in parfiles if x.split('/')[-1].split('.')[0] in PsrList]
timfiles = [x for x in timfiles if x.split('/')[-1].split('_')[0] in PsrList]
noisefiles = [x for x in noisefiles if x.split('/')[-1].split('_')[0] in PsrList]

psrs = []
for p, t in zip(parfiles, timfiles):
    psr = Pulsar(p, t, ephem='DE440')
    psrs.append(psr)
    
# set reference time for the sin wave to the earliest TOA in the data set
dataset_tmin = np.min([p.toas.min() for p in psrs])
dataset_tmax = np.max([p.toas.max() for p in psrs])

# sin wave model
@function
def sine_wave(toas, flags, A = -9, f = -9, phase = 0.0):
    return 10 ** A * np.sin(2 * np.pi * (10 ** f) * (toas - dataset_tmin) + phase)
    
def sine_signal(A, f, phase, selection=Selection(selections.no_selection), name = ""):
    return Deterministic(sine_wave(A = A, f = f, phase = phase), selection = selection, name = name)

if options.sin_wave:
    m1 = sine_signal(A = parameter.Uniform(-9, -4)('common_sin_log10_A'),
                     f = parameter.Uniform(-9, -7.7)('common_sin_log10_f'),
                     phase = parameter.Uniform(0, 2 * np.pi)('common_sin_phase'), name='common_sin')
else:
    m1 = None

# load noise models and files
params = {}
for nf in noisefiles:
    with open(nf, 'r') as fin:
        params.update(json.load(fin))

if not options.red_components:
    try:
        red_dict = {}
        with open(noisedir + '/red_dict.json','r') as rd:
            red_dict.update(json.load(rd))
    except:
        raise UserWarning('Custom pulsar red noise frequency components not set.')
else:
    red_dict = options.red_components

if not options.dm_components:
    try:
        dm_dict = {}
        with open(noisedir + '/dm_dict.json','r') as dd:
            dm_dict.update(json.load(dd))
    except:
        raise UserWarning('Custom pulsar DM noise frequency components not set.')
else:
    dm_dict = options.dm_components

if not options.chrom_components:
    try:
        chrom_dict = {}
        with open(noisedir + '/chrom_dict.json','r') as cd:
            chrom_dict.update(json.load(cd))
    except:
        raise UserWarning('Custom pulsar scattering noise frequency components not set.')
else:
    chrom_dict = options.chrom_components

# setup model
pta = models.model_general(psrs, noisedict=params, orf=options.orf,
                           common_psd=options.common_psd,
                           common_components=options.common_components,
                           bayesephem=options.bayesephem,
                           sat_orb_elements=True, tnequad=True,
                           tm_svd=True, tm_marg=False,
                           red_var=True, red_components=red_dict,
                           dm_var=True, dm_components=dm_dict,
                           chrom_var=True, chrom_components=chrom_dict,
                           chrom_kernel='diag', tndm=True,
                           num_dmdips=options.num_dmdips,
                           dmpsr_list=['J1713+0747'], dm_expdip_idx=[1,4],
                           dm_expdip_tmin=[57490,54650],
                           dm_expdip_tmax=[57530,54850], extra_sigs=m1)

# load CRN chains
chain = np.genfromtxt(outdir + 'chain_1.txt')
names = np.loadtxt(outdir + 'pars.txt', dtype=str)
burn = int(0.25 * chain.shape[0])
chain = chain[burn:]

# run and save OS
orf = ['monopole','dipole','hd']
try:
    gamma_common = float(options.gamma_common)
    outdir += 'os_fixed/'
except:
    gamma_common = None
    outdir += 'os_varied/'
for i in orf:
    Op = op.OptimalStatistic(psrs=psrs, pta=pta, orf=i, gamma_common=gamma_common)
    out = Op.compute_noise_marginalized_os(chain=chain, N=int(options.num))
    np.save(outdir + 'os_{0}.npy'.format(i), out)

Op = op.OptimalStatistic(psrs=psrs, pta=pta, gamma_common=gamma_common)
out = Op.compute_noise_marginalized_multiple_corr_os(chain=chain, N=int(options.num))
np.save(outdir + 'os_mult.npy', out)
