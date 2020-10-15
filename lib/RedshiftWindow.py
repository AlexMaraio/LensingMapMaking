from camb.sources import GaussianSourceWindow


class RedshiftWindow:

    def __init__(self, redshift, sigma, window_type=GaussianSourceWindow):
        self.redshift = redshift
        self.sigma = sigma
        self.source_type = 'lensing'
        self.window_type = window_type

    def construct_camb_instance(self):
        return self.window_type(redshift=self.redshift, source_type=self.source_type, sigma=self.sigma)
