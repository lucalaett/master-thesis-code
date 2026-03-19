import os, sys, json
import numpy

psr_list = numpy.loadtxt('psrlist.txt',dtype=str)
noisedir = 'noisefiles'

try:
    os.makedirs(noisedir)
except:
    pass

red_dict = {}
dm_dict = {}
chrom_dict = {}

for psr_name in psr_list:

    noisedict = {}
    dm_dict.update({"{0}".format(psr_name):None})
    red_dict.update({"{0}".format(psr_name):None})
    chrom_dict.update({"{0}".format(psr_name):None})
    
    with open('{0}-TN.par'.format(psr_name)) as noisefile:
        for l in noisefile:
            if 'TNEF' in l:
                ef = l.split(' ')
                ef[3] = ef[3].replace('\n','')
                noisedict.update({"{0}_{1}_efac".format(psr_name,ef[2]):float(ef[3])})
            elif 'TNEQ' in l:
                eq = l.split(' ')
                eq[3] = eq[3].replace('\n','')
                noisedict.update({"{0}_{1}_log10_tnequad".format(psr_name,eq[2]):float(eq[3])})
            elif 'TNDMGam' in l:
                x = l.split(' ')
                xout = float(x[1].replace('\n',''))
                noisedict.update({"{0}_dm_gp_gamma".format(psr_name):xout})
            elif 'TNDMAmp' in l:
                x = l.split(' ')
                xout = float(x[1].replace('\n',''))
                noisedict.update({"{0}_dm_gp_log10_A".format(psr_name):xout})
            elif 'TNRedGam' in l:
                x = l.split(' ')
                xout = float(x[1].replace('\n',''))
                noisedict.update({"{0}_red_noise_gamma".format(psr_name):xout})
            elif 'TNRedAmp' in l:
                x = l.split(' ')
                xout = float(x[1].replace('\n',''))
                noisedict.update({"{0}_red_noise_log10_A".format(psr_name):xout})
            elif 'TNDMC' in l:
                x = l.split(' ')
                xout = float(x[1].replace('\n',''))
                dm_dict.update({"{0}".format(psr_name):int(xout)})
            elif 'TNRedC' in l:
                x = l.split(' ')
                xout = float(x[1].replace('\n',''))
                red_dict.update({"{0}".format(psr_name):int(xout)})
            elif 'TNScatC' in l:
                x = l.split(' ')
                xout = float(x[1].replace('\n',''))
                chrom_dict.update({"{0}".format(psr_name):int(xout)})
            else:
                pass

    with open('{0}/{1}_noise.json'.format(noisedir,psr_name), 'w') as fout:
        json.dump(noisedict, fout, sort_keys=True, indent=4, separators=(',', ': '))

with open('{0}/dm_dict.json'.format(noisedir), 'w') as fout:
        json.dump(dm_dict, fout, sort_keys=True, indent=4, separators=(',', ': '))
with open('{0}/red_dict.json'.format(noisedir), 'w') as fout:
        json.dump(red_dict, fout, sort_keys=True, indent=4, separators=(',', ': '))
with open('{0}/chrom_dict.json'.format(noisedir), 'w') as fout:
        json.dump(chrom_dict, fout, sort_keys=True, indent=4, separators=(',', ': '))
