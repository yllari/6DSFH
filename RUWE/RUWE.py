import numpy as np
import vaex
import matplotlib.pyplot as plt

def ruwe_separate(df):
    return df[df['ruwe']<1.4], df[df['ruwe']>=1.4]

df = vaex.open("../../6D_full_halo_sample_2_5kpc_no_quality_cuts_wAG.hdf5")
df = df[np.abs(df['MG']) < 1e3]
df = df[np.abs(df['COL_bp_rp']) < 1e3]
df1, df2 = ruwe_separate(df)
max_valMG = np.amax(df2['MG'].values)
min_valMG = np.amin(df2['MG'].values)

max_valCOL = np.amax(df2['COL_bp_rp'].values)
min_valCOL = np.amin(df2['COL_bp_rp'].values)

fig, ax = plt.subplots(1,2)
fig.set_size_inches(8,5)
ax[0].hist2d(df1['COL_bp_rp'].values, df1['MG'].values, bins=100)
ax[0].set_xlim([min_valCOL, max_valCOL])
ax[0].set_ylim([max_valMG, min_valMG])
ax[0].set_xlabel("COL BP-RP")
ax[0].set_ylabel(r"$M_G$")
ax[0].set_title("RUWE<1.4")
ax[0].annotate("N=" + str(df1.count()),xy = (0.1,0.1), color="white")
plt.close()
ax[1].hist2d(df2['COL_bp_rp'].values, df2['MG'].values, bins=100)
ax[1].set_xlim([min_valCOL, max_valCOL])
ax[1].set_ylim([max_valMG, min_valMG])
ax[1].set_xlabel("COL BP-RP")
ax[1].set_title("RUWE>=1.4")
ax[1].annotate("N=" + str(df2.count()),xy = (0.1,0.1), color="white")
fig.suptitle("All Halo dat")
fig.savefig("plot_ruwe.png", dpi=300)
