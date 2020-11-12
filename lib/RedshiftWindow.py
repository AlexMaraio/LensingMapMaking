from camb.sources import GaussianSourceWindow


class RedshiftWindow:

    def __init__(self, redshift, sigma, enable_counts=False, window_type=GaussianSourceWindow):
        self.redshift = redshift
        self.sigma = sigma
        self.window_type = window_type
        self.gal_counts = enable_counts

    def construct_camb_instance(self):
        # If we have galaxy counts to compute, then return two types of source window functions at the same redshift
        if self.gal_counts:
            return [self.window_type(redshift=self.redshift, source_type='counts', sigma=self.sigma),
                    self.window_type(redshift=self.redshift, source_type='lensing', sigma=self.sigma)]

        else:
            return [self.window_type(redshift=self.redshift, source_type='lensing', sigma=self.sigma)]
