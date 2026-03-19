import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import os

def weightedavg(rho, sig):
    weights, avg = 0., 0.
    for r,s in zip(rho, sig):
        weights += 1./(s*s)
        avg += r/(s*s)
        
    return avg/weights, np.sqrt(1./weights)

def calc_HDcc(xi_sorted, rho_sorted, sig_sorted):
    
    npairs = 43
    xi_mean = []
    xi_err = []

    rho_avg = []
    sig_avg = []

    i = 0
    while i < len(xi_sorted):

        xi_mean.append(np.mean(xi_sorted[i:npairs+i]))
        xi_err.append(np.std(xi_sorted[i:npairs+i]))

        r, s = weightedavg(rho_sorted[i:npairs+i], sig_sorted[i:npairs+i])
        rho_avg.append(r)
        sig_avg.append(s)

        i += npairs

    xi_mean = np.array(xi_mean)
    xi_err = np.array(xi_err)
    
    return xi_mean, xi_err, rho_avg, sig_avg

def calc_HDcc_bin(xi_sorted, rho_sorted, sig_sorted):

    xi_bin = np.array([1e-3, 30.0, 50.0, 80.0, 100.0, 120.0, 150.0, 180.0]) * np.pi/180.0
    xi_mean = []
    xi_err = []

    rho_avg = []
    sig_avg = []

    j = 0
    for i in range(len(xi_bin)-1):

        npairs = [x for x in xi_sorted if xi_bin[i] < x < xi_bin[i+1]]

        xi_mean.append(np.mean(npairs))
        npairs = len(npairs)

        xi_err.append(np.std(xi_sorted[j:npairs+j]))

        r, s = weightedavg(rho_sorted[j:npairs+j], sig_sorted[j:npairs+j])
        rho_avg.append(r)
        sig_avg.append(s)

        j += npairs

    xi_mean = np.array(xi_mean)
    xi_err = np.array(xi_err)
    
    return xi_mean, xi_err, rho_avg, sig_avg

def calc_HDcc_2d(xi_sorted, rho_sorted, sig_sorted):
    
    npairs = 43
    xi_mean = []
    xi_err = []

    rho_avg = []
    sig_avg = []

    i = 0
    while i < len(xi_sorted):

        xi_mean.append(np.mean(xi_sorted[i:npairs+i]))
        xi_err.append(np.std(xi_sorted[i:npairs+i]))

        r, s = weightedavg(rho_sorted[:,i:npairs+i].flatten(), sig_sorted[:,i:npairs+i].flatten())
        rho_avg.append(r)
        sig_avg.append(s)

        i += npairs

    xi_mean = np.array(xi_mean)
    xi_err = np.array(xi_err)
    
    return xi_mean, xi_err, rho_avg, sig_avg



datadir = './'
outdir = datadir + 'report/gwb_chain/'
plotdir = datadir + 'plots/os/'

try:
    os.makedirs(plotdir)
except:
    pass

try:
    data_mult = np.load(outdir+'os_mult.npy',allow_pickle=True)
    amp_mult = []
    sn_mult = []
    for i in range(3):
        amp_mult.append(data_mult[3][:,i])
        sn_mult.append(data_mult[3][:,i]/data_mult[4][:,i])
except:
    pass

# load OS output files

amp = []
sn = []

data = np.load(outdir+'os_monopole.npy',allow_pickle=True)
amp.append(data[3])
sn.append(data[4])

data = np.load(outdir+'os_dipole.npy',allow_pickle=True)
amp.append(data[3])
sn.append(data[4])

data = np.load(outdir+'os_hd.npy',allow_pickle=True)
amp.append(data[3])
sn.append(data[4])

idx_sort = np.argsort(data[0])
phi_sort = data[0][idx_sort]
av_sort = data[1][:,idx_sort]
sig_sort = data[2][:,idx_sort]

av_marg = []
av_marg_unc = []
for i in range(len(av_sort)):
    phi_plot, phi_plot_unc, av_plot, av_plot_unc = calc_HDcc_bin(phi_sort,av_sort[i],sig_sort[i])
    av_marg.append(av_plot)
    av_marg_unc.append(av_plot_unc)
av_marg = np.atleast_2d(av_marg)
av_marg_unc = np.atleast_2d(av_marg_unc)

omc2 = (1 - np.cos(phi_plot)) / 2
hd = 1.5 * omc2 * np.log(omc2) - 0.25 * omc2 + 0.5

av_out = []
av_out_unc = []
chi2 = []
bin_size = len(phi_plot)
for i in range(bin_size):
    av_t, av_t_unc = scipy.stats.distributions.norm.fit(av_marg[:,i]/np.median(amp[2]))
    av_out.append(av_t)
    av_out_unc.append(av_t_unc)
    chi2.append((av_t - hd[i])**2./av_t_unc**2.)

# save chi2 distribution
np.savetxt(plotdir+'mean_dev_{0}.txt'.format(bin_size),np.vstack((phi_plot,av_out,av_out_unc,chi2)).T)

# plot OS correlation figure
fig = plt.figure()

plt.errorbar(phi_plot,av_out,xerr=phi_plot_unc,yerr=av_out_unc,capsize=10)

plt.axhline(np.median(amp[0])/np.median(amp[2]),c='k',lw=2, label='monopole')

omc2 = (1 - np.cos(phi_sort)) / 2
hd = 1.5 * omc2 * np.log(omc2) - 0.25 * omc2 + 0.5
plt.plot(phi_sort,hd,'r',lw=2, label='HD')

plt.legend(loc=0)
plt.grid()
plt.xlabel(r'$\theta$ [rad]',fontsize='x-large')

fig.savefig(plotdir+'os_unc_{0}.pdf'.format(bin_size),bbox_inches='tight')
fig.savefig(plotdir+'os_unc_{0}.png'.format(bin_size),bbox_inches='tight')
plt.clf()



orf = ['monopole','dipole','hd']
Colors = ['C0','C1','C2']

# plot amplitude
fig = plt.figure()
r = (-0.5e-29,2.0e-29)

med_amp = []
for i in range(3):
    med,low,high = np.percentile(amp[i],(50,5,95))
    med_amp.append(med)
    title = r"${{{0:.2e}}}_{{-{1:.1e}}}^{{+{2:.1e}}}$"
    title = title.format(med,med-low,high-med)

    plt.hist(amp[i],25,range=r,density=True,histtype='step',lw=2,color=Colors[i],label='{0}: {1}'.format(orf[i],title))

    try:
        med,low,high = np.percentile(amp_mult[i],(50,5,95))
        title = r"${{{0:.2e}}}_{{-{1:.1e}}}^{{+{2:.1e}}}$"
        title = title.format(med,med-low,high-med)
        plt.hist(amp_mult[i],25,range=r,density=True,histtype='step',ls='--',lw=2,color=Colors[i],label='{0}: {1}'.format(orf[i],title))
    except:
        pass

plt.legend(loc=0,fontsize='x-large')
plt.xlabel(r'$A_{CRS}^2$',fontsize='xx-large')
plt.tick_params(labelsize='x-large')

fig.savefig(plotdir+'os_amp.pdf',bbox_inches='tight')
fig.savefig(plotdir+'os_amp.png',bbox_inches='tight')
plt.clf()

# plot S/N
fig = plt.figure()
r = (-2,8)

med_sn = []
for i in range(3):
    med,low,high = np.percentile(sn[i],(50,5,95))
    med_sn.append(med)
    title = r"${{{0:.2f}}}_{{-{1:.1f}}}^{{+{2:.1f}}}$"
    title = title.format(med,med-low,high-med)

    plt.hist(sn[i],30,range=r,density=True,histtype='step',lw=2,color=Colors[i],label='{0}: {1}'.format(orf[i],title))

    try:
        med,low,high = np.percentile(sn_mult[i],(50,5,95))
        title = r"${{{0:.2f}}}_{{-{1:.1f}}}^{{+{2:.1f}}}$"
        title = title.format(med,med-low,high-med)
        plt.hist(sn_mult[i],30,range=r,density=True,histtype='step',ls='--',lw=2,color=Colors[i],label='{0}: {1}'.format(orf[i],title))
    except:
        pass

plt.legend(loc=0,fontsize='x-large')
plt.xlabel('S/N',fontsize='xx-large')
plt.tick_params(labelsize='x-large')

fig.savefig(plotdir+'os_sn.pdf',bbox_inches='tight')
fig.savefig(plotdir+'os_sn.png',bbox_inches='tight')
plt.clf()
