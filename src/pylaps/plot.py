import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def plot_1d_cut(xgrid, ygrid, data, time, ind, axis, direction, nt):
    """
    Plot 1D cut along a specified direction (x or y) and save the figure.

    Parameters:
    xgrid, ygrid (np.ndarray): Grid points.
    data (dict): Dictionary of data arrays.
    time (float): Time value.
    ind (int): Index for slicing the data.
    axis (str): Axis along which to plot ('x' or 'y').
    direction (str): Direction of the cut ('along_x' or 'along_y').
    nt (int): Time step number.
    """
    fig = plt.figure(figsize=[10, 10])

    for i, (key, values) in enumerate(data.items()):
        sub = fig.add_subplot(4, 1, i+1)
        if axis == 'x':
            sub.plot(xgrid, values[:, ind], label=key)
        else:
            sub.plot(ygrid, values[ind, :], label=key)
        sub.legend()

    fig.suptitle(r'$t={:.2f}$'.format(time))
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f'./figure/cut_{direction}_{nt:03d}.png')
    plt.close(fig)

def plot_2d_map(XX, YY, data, title, vlim=None, if_div_cmap=True):
    """
    Plot 2D color-coded maps and save the figure.

    Parameters:
    XX, YY (np.ndarray): Meshgrid arrays.
    data (np.ndarray): Data array to be plotted.
    title (str): Title of the plot.
    filename (str): Filename to save the plot.
    vlim (list, optional): Value limits for the color map.
    if_div_cmap (bool, optional): Whether to use a diverging color map.
    """
    fig = plt.figure()
    sub = fig.add_subplot(111)
    if if_div_cmap:
        pm = sub.pcolormesh(XX, YY, data, shading='gouraud', cmap='RdBu', norm=TwoSlopeNorm(vcenter=0))
    else:
        if vlim is None:
            pm = sub.pcolormesh(XX, YY, data, shading='gouraud', cmap=plt.cm.plasma)
        else:
            pm = sub.pcolormesh(XX, YY, data, shading='gouraud', cmap=plt.cm.plasma, vmin=vlim[0], vmax=vlim[1])
    
    sub.set_aspect('equal')
    fig.colorbar(pm, shrink=0.9, aspect=15)
    sub.set_xlabel(r'$x$', fontsize=12)
    sub.set_ylabel(r'$y$', fontsize=12)
    sub.set_title(title)
    return fig, sub