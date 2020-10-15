import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=1.75, rc={'text.usetex': True})

from lib.RedshiftWindow import RedshiftWindow
from lib.CambObject import CambObject

import camb

if __name__ == '__main__':
    window1 = RedshiftWindow(0.2, 0.05)
    window2 = RedshiftWindow(0.49, 0.05)
    window3 = RedshiftWindow(0.62, 0.05)
    window4 = RedshiftWindow(0.74, 0.05)
    window5 = RedshiftWindow(0.85, 0.05)
    window6 = RedshiftWindow(0.96, 0.05)
    window7 = RedshiftWindow(1.08, 0.05)
    window8 = RedshiftWindow(1.24, 0.05)
    window9 = RedshiftWindow(1.45, 0.05)
    window10 = RedshiftWindow(2.04, 0.15)

    window1 = RedshiftWindow(0.5, 0.05)
    window2 = RedshiftWindow(2.0, 0.05)

    linear = CambObject('Linear', 5000, non_linear=False)
    non_linear = CambObject('Non-linear', 5000, non_linear=True)

    linear.set_window_functions([window1, window2, window10])
    # non_linear.set_window_functions([window1, window5, window10])
    non_linear.set_window_functions([window1, window2])

    c_ells_linear = linear.get_c_ells_dict()
    c_ells_nonlinear = non_linear.get_c_ells_dict()

    plt.figure(figsize=(13, 7))
    plt.loglog(non_linear.ells, c_ells_nonlinear['W1xW1'][2:], color='r', label=r'$1 \times 1$')
    plt.loglog(non_linear.ells, c_ells_nonlinear['W2xW2'][2:], color='y', label=r'$2 \times 2$')
    # plt.loglog(non_linear.ells, c_ells_nonlinear['W3xW3'][2:], color='b', label=r'$10 \times 10$')
    plt.loglog(linear.ells, c_ells_linear['W1xW1'][2:], color='r', ls='--')
    plt.loglog(linear.ells, c_ells_linear['W2xW2'][2:], color='y', ls='--')
    # plt.loglog(linear.ells, c_ells_linear['W3xW3'][2:], color='b', ls='--')
    plt.title('Lensing power spectrum')
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
    plt.legend()
    plt.tight_layout()
    plt.show(block=False)

    non_linear.output_c_ells()
    non_linear.split_camb_output()

    # non_linear.plot_1x2_map_power_spectrum(key='W1xW1', nside=2048)

    non_linear.create_fields_info()

    non_linear.write_flask_config_file()

    non_linear.set_flask_executable('~/Documents/PhD/Codes/flask/bin/flask')

    non_linear.run_flask()

    # ! non_linear.estimate_cl_from_alm()

    non_linear.trim_flask_alm_output()
    non_linear.use_cpp_thingy()

    non_linear.plot_flask_output()
