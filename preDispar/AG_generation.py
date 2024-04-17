### FUNCTION!
# Reddening computation:
def reddening_computation(self, dust_map):

    # Compute ebv from one of the published 3D maps:

    if dust_map == 'l22':
        # We use the dust map from Lallement 2022.

        EBV = np.array(l22.A0(self.evaluate('l'),self.evaluate('b'),self.evaluate('distance')*1000.))/3.1

        l = np.asarray(self.evaluate('l'))
        b = np.asarray(self.evaluate('b'))
        r = np.asarray(self.evaluate('distance')*1000.)

        lrad=np.radians(l)
        brad=np.radians(b)

        x=r*np.cos(brad)*np.cos(lrad)
        y=r*np.cos(brad)*np.sin(lrad)
        z=r*np.sin(brad)

        err_distance=self.evaluate('distance')/self.evaluate('parallax_over_error')
        err_distance[~np.isfinite(err_distance)] = np.max(err_distance[np.isfinite(err_distance)])

        iout=np.where((np.abs(x) >= 5000) | (np.abs(y) >= 5000) | (np.abs(z) >= 400))

        A0sup=l22.A0(l,b,self.evaluate('distance')*1000.+err_distance*1000.)
        A0inf=l22.A0(l,b,self.evaluate('distance')*1000.-err_distance*1000.)
        delta_A0=np.sqrt( (3.1*EBV*0.01)**2 + (np.abs(np.array(A0sup)-np.array(A0inf))/2)**2 )

        if x[iout].any(): delta_A0[iout] = 3.1*EBV[iout]*0.01  ## for a star beyond the limits of the extinction map we consider that the error in the extinction is only affected by the error in the methodology to infer A0 (a few percents ~ 0.01*A0). For the other cases the error is dominated by the error in the distance assume for the star which is where we calculate the A0 in the map.

        delta_EBV = delta_A0/3.1

        del iout, lrad, brad, x, y, z, A0sup, A0inf


    elif dust_map == 'l18':
        # We use the dust map from lallement18. It uses dust_maps_3d/lallement18.py and the file lallement18_nside64.fits.gz
        EBV = l18.ebv(self.evaluate('l'), self.evaluate('b'), self.evaluate('distance')*1000.0)
        delta_EBV = EBV*0.0

    elif dust_map == 'p19':
        # We use the dust map from pepita19. It would use dust_maps_3d/pepita19.py and the file pepita19.fits.gz
        EBV = p19.ebv(self.evaluate('l'), self.evaluate('b'), self.evaluate('distance'))
        delta_EBV = EBV*0.0

    elif dust_map == 'bayestar':

        err_distance=self.evaluate('distance')/self.evaluate('parallax_over_error')
        err_distance[~np.isfinite(err_distance)] = np.max(err_distance[np.isfinite(err_distance)])

        sup_dist, inf_dist = self.evaluate('distance')+err_distance, self.evaluate('distance')-err_distance

        bayestar = BayestarQuery(max_samples=5, version='bayestar2019')
        coords = SkyCoord(self.evaluate('l')*units.deg, self.evaluate('b')*units.deg, distance=self.evaluate('distance')*units.kpc, frame='galactic')
        coordssup = SkyCoord(self.evaluate('l')*units.deg, self.evaluate('b')*units.deg, distance=abs(sup_dist*units.kpc), frame='galactic') # In the simulation of distances, some might go below 0, I make them positive.
        coordsinf = SkyCoord(self.evaluate('l')*units.deg, self.evaluate('b')*units.deg, distance=abs(inf_dist*units.kpc), frame='galactic') # In the simulation of distances, some might go below 0, I make them positive.

        EBV = bayestar(coords, mode='mean')
        EBVsup=bayestar(coordssup, mode='mean')
        EBVinf=bayestar(coordsinf, mode='mean')

        EBV_perc = bayestar(coords, mode='percentile', pct=[16., 84.])
        delta_EBV = np.sqrt(((EBV_perc[:,1]-EBV_perc[:,0])/2.0)**2 + (np.abs(np.array(EBVsup)-np.array(EBVinf))/2)**2)

        del bayestar, EBV_perc

        # Problems with Green catalogue: 1- Not EBV; 2- zero values.
        # 1: True_EBV = 0.884*EBV, taken from https://urldefense.com/v3/__http://argonaut.skymaps.info/usage*units__;Iw!!D9dNQwwGXtA!XAuX6itZfdSRbbCZALx5jRVK7sEXqqMv1MCsMGaJDsR6eYEYqoGACWehcWclRLc8mmaQLsIk4liQun095lxYTsSptQTT$
        EBV, delta_EBV = 0.884*EBV, 0.884*delta_EBV
        # 2: Discretization from Green email, if EBV = 0.0, assume a number in the range 0.0-0.02 (positive I guess, right?)
        #coords = np.where(EBV==0.0)
        #EBV[coords] = np.random.uniform(0.0, 0.02, len(EBV[coords]))
        #coords = np.where(delta_EBV==0.0)
        #delta_EBV[coords] = np.random.uniform(0.0, 0.02, len(delta_EBV[coords]))

    # From those E(B-V) and errors, compute the AG and E_BP_RP values (several recipes in the literature)

    if ext_rec == 'f19':
        ### Fitzpatrick + 2019 ---> https://urldefense.com/v3/__https://www.cosmos.esa.int/web/gaia/edr3-extinction-law__;!!D9dNQwwGXtA!XAuX6itZfdSRbbCZALx5jRVK7sEXqqMv1MCsMGaJDsR6eYEYqoGACWehcWclRLc8mmaQLsIk4liQun095lxYTuOQaWja$

        A0 = 3.1*EBV
        bp_rp_0 = self.evaluate('bp_rp')
        bp_rp = self.evaluate('bp_rp')
        # We enter in an iterative process of 10 steps to refine the kBP, kRP, and kG values. 'Cause we do not have intrinsic colour.
        for i in range(5):
            kG  =( 0.996727636184546 -0.161532912039159*bp_rp_0 + 0.0121894634053956*bp_rp_0**2 + 0.00108655443373807*bp_rp_0**3
                 -0.0378418720324357*A0 + 0.00150218619071414*A0**2
                 -2.43558666819941e-05*A0**3 + 0.0119452762544019*bp_rp_0*A0
                 -0.000900410348061063*A0*bp_rp_0**2
                 -0.000285148019567325*bp_rp_0*A0**2)

            kBP =( 1.15477509956089 -0.0860151175561996*bp_rp_0
                  -0.0302071198654238*bp_rp_0**2 + 0.0140420820336445*bp_rp_0**3
                  -0.0223968266786066*A0 + 0.000823544766886717*A0**2
                  -1.20366204386921e-05*A0**3 + 0.00699684361199003*bp_rp_0*A0
                  -0.000499920649896623*A0*bp_rp_0**2
                  -0.000148686538559682*bp_rp_0*A0**2)

            kRP = (0.662595389120923 -0.0167121614962357*bp_rp_0
                 - 0.00144198404916194*bp_rp_0**2 -0.00133749382473256*bp_rp_0**3
                 - 0.00655783097133718*A0 + 3.36125225635362e-05*A0**2
                 + 1.61356747065676e-06*A0**3 + 2.20748538619438e-05*bp_rp_0*A0
                 + 0.000178453487002225*A0*bp_rp_0**2 + 1.0336186431516e-05*bp_rp_0*A0**2)

            bp_rp_0 = bp_rp - (kBP-kRP)*A0 # 0 -- intrinsic, otherwise observed, Ax = mx - m0x --> m0x = mx - Ax = mx - kx*A0  ==> m0BP - m0RP = mbp - mrp - (kbp - krp)*A0 = bp_rp

        AG = kG*A0
        Abp = kBP*A0
        Arp = kRP*A0
        E_BP_RP = (kBP-kRP)*A0 # A_BP - A_RP

        delta_A0 = 3.1*delta_EBV
        delta_kG = (-0.0378418720324357 + 2*0.00150218619071414*A0
                    -3*2.43558666819941e-05*A0**2 + 0.0119452762544019*bp_rp_0
                    -0.000900410348061063*bp_rp_0**2 - 2*0.000285148019567325*bp_rp_0*A0)*delta_A0  #### derivada con respect A0 * delta_a0
        delta_AG = np.sqrt((A0*delta_kG)**2 + (kG*delta_A0)**2)

        delta_kBP_kRP = (- 0.015838995707269418 + 0.0015798644886463617*A0 - 4.0950563728046576e-05*A0**2
                         + 0.006974768758128086*bp_rp_0 - 0.000678374136898848*bp_rp_0**2
                         - 0.000318045449982396*bp_rp_0*A0 )*delta_A0  ##### restar KBP-kRP y derivar respecto a A0
        delta_E_BP_RP = np.sqrt((A0*delta_kBP_kRP)**2 + ((kBP-kRP)*delta_A0)**2)
        del A0, bp_rp_0, kG, kBP, kRP, delta_kBP_kRP, delta_A0, delta_kG

    elif ext_rec == 'cas18':

        # Colour--Teff relation determined from Gaia DR2 data --> Computed using color_teff_GDR2.ipynb
        poly = np.array([62.55257114, -975.10845442, 4502.86260828, -8808.88745592, 10494.72444183])
        Teff = np.poly1d(poly)(self.evaluate('bp_rp')) / 1e4

        # We enter in an iterative process of 10 steps to refine the AG and E_BP_RP computation. 'Cause we do not have intrinsic colour.
        for i in range(10):
            E_BP_RP = (-0.0698 + Teff*(3.837-2.530*Teff)) * EBV # A_BP - A_RP
            AG = (1.4013 + Teff*(3.1406-1.5626*Teff)) * EBV
            Teff = np.poly1d(poly)(self.evaluate('bp_rp') - E_BP_RP) / 1e4
        delta_AG = (1.4013 + Teff*(3.1406-1.5626*Teff)) * delta_EBV
        delta_E_BP_RP = (-0.0698 + Teff*(3.837-2.530*Teff)) * delta_EBV

    elif ext_rec == 'bab18':
        ### Babusiaux et al. 2018 --> https://urldefense.com/v3/__http://aselfabs.harvard.edu/abs/2018A*26A...616A..10G__;JQ!!D9dNQwwGXtA!XAuX6itZfdSRbbCZALx5jRVK7sEXqqMv1MCsMGaJDsR6eYEYqoGACWehcWclRLc8mmaQLsIk4liQun095lxYTnJ0NW5e$  (eq 1, Table 1)

        A0 = 3.1*EBV
        bp_rp_0 = self.evaluate('bp_rp')
        bp_rp = self.evaluate('bp_rp')

        # We enter in an iterative process of 10 steps to refine the kBP, kRP, and kG values. 'Cause we do not have intrinsic colour.
        for i in range(3):
            kG  = 0.9761 - 0.1704*bp_rp_0 + 0.0086*bp_rp_0**2 + 0.0011*bp_rp_0**3 - 0.0438*A0 + 0.00130*A0**2 + 0.0099*bp_rp_0*A0
            kBP = 1.1517 - 0.0871*bp_rp_0 - 0.0333*bp_rp_0**2 + 0.0173*bp_rp_0**3 - 0.0230*A0 + 0.00060*A0**2 + 0.0043*bp_rp_0*A0
            kRP = 0.6104 - 0.0170*bp_rp_0 - 0.0026*bp_rp_0**2 - 0.0017*bp_rp_0**3 - 0.0078*A0 + 0.00005*A0**2 + 0.0006*bp_rp_0*A0
            bp_rp_0 = bp_rp_0 - (kBP-kRP)*A0 # 0 -- intrinsic, otherwise observed, Ax = mx - m0x --> m0x = mx - Ax = mx - kx*A0  ==> m0BP - m0RP = mbp - mrp - (kbp - krp)*A0 = bp_rp

        AG = kG*A0
        delta_A0 = 3.1*delta_EBV
        delta_kG = (-0.0438 + 2.0*0.00130*A0 + 0.0099*bp_rp_0)*delta_A0
        delta_AG = np.sqrt((A0*delta_kG)**2 + (kG*delta_A0)**2)
        del delta_kG

        E_BP_RP = (kBP-kRP)*A0 # A_BP - A_RP
        delta_kBP_kRP = (-0.0152 + 2.0*0.00055*A0 + 0.0037*bp_rp_0)*delta_A0
        delta_E_BP_RP = np.sqrt((A0*delta_kBP_kRP)**2 + ((kBP-kRP)*delta_A0)**2)
        del A0, bp_rp_0, kG, kBP, kRP, delta_kBP_kRP, delta_A0

    elif ext_rec == 'else':
        AG = 0.86117*3.1*EBV
        Abp = 1.06126*3.1*EBV
        Arp = 0.64753*3.1*EBV
        E_BP_RP = Abp - Arp



    COL_bp_rp = self.evaluate('bp_rp') - E_BP_RP
    MG = self.evaluate('phot_g_mean_mag')-AG+5.-5.*np.log10(1000.0/self.evaluate('new_parallax'))

    # Cleaning values
    thresh = 100
    E_BP_RP[E_BP_RP > thresh] = np.nan
    delta_E_BP_RP[delta_E_BP_RP > thresh] = np.nan
    delta_AG[delta_AG > thresh] = np.nan
    COL_bp_rp[COL_bp_rp > thresh] = np.nan
    MG[MG > thresh] = np.nan
    return E_BP_RP, delta_E_BP_RP, AG, delta_AG, COL_bp_rp, MG

## PREAMBLE OF CODE:
import astropy.units as units
from astropy.coordinates import SkyCoord
import L22Map as l22
#import lallement18_distant_all_gaia as l18
from dustmaps.config import config
config.reset()
config['data_dir'] = '/scratch1/truiz/cylindre_nov23/selections_and_catalogues/dustmaps'
import dustmaps.bayestar
dustmaps.bayestar.fetch()
from dustmaps.bayestar import BayestarQuery

dust_map = 'bayestar' # bayestar, l22
ext_rec = 'f19'

## Adding this info to an hdf5 file:
gaia = vaex.open('6D_full_halo_sample_2_5kpc_no_quality_cuts.hdf5')
E_BP_RP2, delta_E_BP_RP2, AG2, delta_AG2, COL_bp_rp2, MG2 = reddening_computation(gaia, dust_map)
gaia.add_column('E_BP_RP',E_BP_RP2)
gaia.add_column('delta_E_BP_RP',delta_E_BP_RP2)
gaia.add_column('AG',AG2)
gaia.add_column('delta_AG',delta_AG2)
gaia.add_column('COL_bp_rp',COL_bp_rp2)
gaia.add_column('MG',MG2)
gaia.export_hdf5("6D_full_halo_sample_2_5kpc_no_quality_cuts_wAG.hdf5", progress=True)
del E_BP_RP2, delta_E_BP_RP2, AG2, delta_AG2, COL_bp_rp2, MG2
