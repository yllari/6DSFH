import vaex
import numpy as np
class PreDisPar():
    """Data procesing class in preparation for disPar

    When instantiated, a file to read, group and columns (optional) are given
    """
    data = None
    sel_group = None

    def __init__(self, in_file, group, columns=None):
        self.data = vaex.open(in_file)
        self.sel_group = group
        self.data = self.data[self.data['derived_labels_group'] == group]

        if columns != None:
            self.data = self.data[columns]

    def get_full(self):
        """Process data as with given critera.

        A function should be defined by critera
        """
        # -+-+-+-+-+- Filters -+-+-+-+-+-+-+-
        # General
        bright = 'MG<5.5'
        dist = '1/new_parallax < 1.2'
        df = self.data.filter(bright).extract()
        df.select(dist, name="dist")
        return df.filter('dist').extract()

    def get_qshag(self):
        """Process data as QSHAG standard

        Filters in extinction, parallax over error to
        ensure dist = 1/parallax and quality with excess color
        """
        # QSHAG
        bright = 'MG<5.5'
        extinction = '(AG < 0.5)'
        poege = '(parallax_over_error > 5.0)'
        quality =  '(0.001+0.039*(bp_rp) < log10(phot_bp_rp_excess_factor)) & (log10(phot_bp_rp_excess_factor) < 0.12 + 0.039*(bp_rp))'
        dist = '1/new_parallax < 1.2'
        vel_cut = 'radial_velocity_error < 20'

        df = self.data.filter(bright).extract()
        # Separating just for the sake of maintaining
        # self.data selections clean
        df.select(extinction, name="ext")
        df.select(poege, name="poege")
        df.select(quality, name="quality")
        df.select(dist, name="dist")
        df.select(vel_cut, name="vel_cut")
        return df.filter('(ext)&(poege)&(quality)&(dist)&(vel_cut)').extract()
    def print_data(self):
        print(self.data)
        return None


    def print_columns(self):
        print(self.data.column_names)
        return None


if __name__ == "__main__":
    in_file = "../6D_full_halo_sample_2_5kpc_no_quality_cuts_wAG.hdf5"
    g = vaex.open(in_file)
    # Reduced dictionary with mapping of labels to names
    lb_dc= {
        "Smooth": -1,
        "GE": 0,
        "LRL3": 1,
        "HotDsk": 2,
        "Thamnos": 3,
        "Helmi": 4,
        "Sequoia": 5,
    }
    group = lb_dc["GE"]


    pre_dispar = PreDisPar(in_file, group)
    out_data = pre_dispar.get_qshag()
    print("number of stars", out_data.count())
    out_data.export_hdf5('GE_qshag.hdf5', progress=True)
