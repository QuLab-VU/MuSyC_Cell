import pandas as pd
from pylab import *
from mpl_toolkits.mplot3d import *

from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm

from synergy_tools import *

def plot_surface(d1_min, d1_max, d2_min, d2_max, E0, E1, E2, E3, h1, h2, r1,   \
                 r1r, r2, r2r, alpha1, alpha2, scatter_points=None, N=50,      \
                 elev=20, azim=19, fname=None, zlim=None, cmap='bwr',          \
                 xlabel="x", ylabel="y", zlabel="z", d1buffer=0., d2buffer=0., \
                 zslice=None):
    
    fig = figure(figsize=(11,6))
    ax = fig.gca(projection='3d')
    
    d1 = logspace(d1_min-d1buffer, d1_max, N)
    d2 = logspace(d2_min-d2buffer, d2_max, N)
    DD1, DD2 = meshgrid(d1, d2)
    E = hill_2D(DD1, DD2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)

    # Used for colorbar
    zmin = np.min(E)
    zmax = np.max(E)
    if zslice is not None:
        top = np.abs(zmax-zslice)
        bot = np.abs(zmin-zslice)
        m = max(top,bot)
        vmin, vmax = zslice-m, zslice+m
    else:
        vmin, vmax = None, None

    # Plot it once on a throwaway axis, to get cbar without alpha problems
    my_cmap_rgb = plt.get_cmap(cmap)(np.arange(256))
    alpha = 0.8

    for i in range(3): # Do not include the last column!
        my_cmap_rgb[:,i] = (1 - alpha) + alpha*my_cmap_rgb[:,i]
    my_cmap = ListedColormap(my_cmap_rgb, name='my_cmap')

    surf = ax.plot_surface(np.log10(DD1), np.log10(DD2), E, cstride=1, rstride=1, cmap = my_cmap, vmin=vmin, vmax=vmax, linewidth=0)

    cbar = fig.colorbar(surf)
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor("face")

    cla()


    # Plot the surface for real, with appropriate alpha
    surf = ax.plot_surface(np.log10(DD1), np.log10(DD2), E, cstride=1, rstride=1, alpha=alpha, cmap = cm.bwr, vmin=vmin, vmax=vmax, linewidth=0)


    
    # colored curves on left and right
    lw = 5
    E1_alone = hill_2D(d1, np.power(10.,d2_min-d2buffer), E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    E2_alone = hill_2D(np.power(10.,d1_min-d1buffer), d2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    
    E1_with_2 = hill_2D(d1, np.power(10.,d2_max), E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    E2_with_1 = hill_2D(np.power(10.,d1_max), d2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
    
    ax.plot(np.log10(d1), d2_min*ones(d1.shape)-d2buffer, E1_alone, linewidth=lw)
    ax.plot(d1_min*ones(d2.shape)-d1buffer, np.log10(d2), E2_alone, linewidth=lw)
    
    ax.plot(np.log10(d1), d2_max*ones(d1.shape), E1_with_2, linewidth=lw)
    ax.plot(d1_max*ones(d2.shape), np.log10(d2), E2_with_1, linewidth=lw)
    
    
    # light grey grid across surface
    for y in linspace(d2_min-d2buffer,d2_max,N/5):
        E_slice = hill_2D(d1, np.power(10.,y), E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
        ax.plot(np.log10(d1), y*np.ones(d2.shape), E_slice, '-k', linewidth=2, alpha=0.2)


    for x in linspace(d1_min-d1buffer, d1_max,N/5):
        E_slice = hill_2D(np.power(10.,x), d2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
        ax.plot(x*np.ones(d1.shape), np.log10(d2), E_slice, '-k', linewidth=2, alpha=0.2)


    # Set the view
    ax.view_init(elev=elev, azim=azim)


    # Label axes appropriately
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)


    if scatter_points is not None and type(scatter_points)==pd.DataFrame:
        s_d1 = np.log10(scatter_points['drug1.conc'])
        s_d2 = np.log10(scatter_points['drug2.conc'])
        s_d1.loc[s_d1==-np.inf] = d1_min-d1buffer
        s_d2.loc[s_d2==-np.inf] = d2_min-d2buffer
        
        ax.scatter(s_d1, s_d2, scatter_points['effect'], depthshade=True)

        # Plot error bars        
        if "effect.95ci" in scatter_points.columns:
            for _d1, _d2, _E, _erb in zip(s_d1, s_d2, scatter_points['effect'], scatter_points['effect.95ci']):
                ax.plot([_d1,_d1], [_d2,_d2], [_E-_erb, _E+_erb], 'k-', alpha=0.3)


    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    if zslice is not None:
        # Plane of DIP=0
        c_plane = colorConverter.to_rgba('k', alpha=0.2)
        verts = [array([(d1_min-d1buffer,d2_min-d2buffer), (d1_min-d1buffer,d2_max), (d1_max,d2_max), (d1_max,d2_min-d2buffer), (d1_min-d1buffer,d2_min-d2buffer)])]
        poly = PolyCollection(verts, facecolors=c_plane)
        ax.add_collection3d(poly, zs=[zslice], zdir='z')

        # Plot intersection of surface with DIP=0 plane
        CS = contour(np.log10(DD1),np.log10(DD2),E,levels=[zslice], linewidths=3, colors='k')

    # If needed, manually set zlim
    if zlim is not None: ax.set_zlim(zlim)

    # Save plot, or show it
    if fname is None: show()
    else: plt.savefig(fname)
    clf()
    plt.close()
    


def plot_arbitrary_surface(DD1, DD2, E, scatter_points=None, N=50, elev=20,   \
                 azim=19, fname=None, zlim=None, cmap='bwr', xlabel="x",      \
                 ylabel="y", zlabel="z", title=None, zslice=None, edges=True, \
                 logspace=False):
    """
    DD1 and DD2 are a meshgrid defining the domain of E
    scatter_points is an optional pandas DataFrame with individual datapoints to scatter
        drug1.conc
        drug2.conc
        effect
        effect.95ci (*optional*)
    zslice is an optional float that will highlight a specific xy-plane along z-axis
    edges is an optional bool. If edges=True, the 4 edges of the surface will be highlighted
    logsace is an optional bool. If logspace=True, DD1 and DD2 are expected to be
        logtransformed concentrations. Otherwise they are assumed to be in linear space,
        and will be transformed to log values for plotting.
    """
    fig = figure(figsize=(3.5,3))
    ax = fig.gca(projection='3d')
    
    
    if logspace:
        d1_min = np.min(DD1)
        d1_max = np.max(DD1)
        d2_min = np.min(DD2)
        d2_max = np.max(DD2)    
        DD1 = np.power(10., DD1)
        DD2 = np.power(10., DD2)

    else:
        d1_min = np.log10(np.min(DD1))
        d1_max = np.log10(np.max(DD1))
        d2_min = np.log10(np.min(DD2))
        d2_max = np.log10(np.max(DD2))
        
    d1 = np.logspace(d1_min, d1_max, DD1.shape[1])
    d2 = np.logspace(d2_min, d2_max, DD2.shape[0])

    # Used for colorbar
    zmin = np.min(E)
    zmax = np.max(E)
    if zslice is not None:
        top = np.abs(zmax-zslice)
        bot = np.abs(zmin-zslice)
        m = max(top,bot)
        vmin, vmax = zslice-m, zslice+m
    else:
        vmin, vmax = None, None

    # Plot it once on a throwaway axis, to get cbar without alpha problems
    my_cmap_rgb = plt.get_cmap(cmap)(np.arange(256))
    alpha = 0.8

    for i in range(3): # Do not include the last column!
        my_cmap_rgb[:,i] = (1 - alpha) + alpha*my_cmap_rgb[:,i]
    my_cmap = ListedColormap(my_cmap_rgb, name='my_cmap')

    surf = ax.plot_surface(np.log10(DD1), np.log10(DD2), E, cstride=1, rstride=1, cmap = my_cmap, vmin=vmin, vmax=vmax, linewidth=0)

    cbar = fig.colorbar(surf)
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor("face")

    cla()


    # Plot the surface for real, with appropriate alpha
    surf = ax.plot_surface(np.log10(DD1), np.log10(DD2), E, cstride=1, rstride=1, alpha=alpha, cmap = cm.bwr, vmin=vmin, vmax=vmax, linewidth=0)

    # colored curves on left and right
    lw = 5
    E1_alone = E[0,:]
    E2_alone = E[:,0]
    
    E1_with_2 = E[-1,:]
    E2_with_1 = E[:,-1]
    
    ax.plot(np.log10(d1), d2_min*ones(d1.shape), E1_alone, linewidth=lw)
    ax.plot(d1_min*ones(d2.shape), np.log10(d2), E2_alone, linewidth=lw)
    
    ax.plot(np.log10(d1), d2_max*ones(d1.shape), E1_with_2, linewidth=lw)
    ax.plot(d1_max*ones(d2.shape), np.log10(d2), E2_with_1, linewidth=lw)

    # Set the view
    ax.view_init(elev=elev, azim=azim)

    # Label axes appropriately
    ax.set_ylabel("log(%s)[M]"%drug1_name,labelpad=-5)   
    ax.set_xlabel("log(%s)[M]"%drug2_name,labelpad=-5)
    ax.set_zlabel(metric_name,labelpad=-3)

    if scatter_points is not None and type(scatter_points)==pd.DataFrame:
        s_d1 = np.log10(scatter_points['drug1.conc'])
        s_d2 = np.log10(scatter_points['drug2.conc'])
        s_d1.loc[s_d1==-np.inf] = d1_min
        s_d2.loc[s_d2==-np.inf] = d2_min
        
        ax.scatter(s_d1, s_d2, scatter_points['effect'], depthshade=True)

        # Plot error bars        
        if "effect.95ci" in scatter_points.columns:
            for _d1, _d2, _E, _erb in zip(s_d1, s_d2, scatter_points['effect'], scatter_points['effect.95ci']):
                ax.plot([_d1,_d1], [_d2,_d2], [_E-_erb, _E+_erb], 'k-', alpha=0.3)


    #format axis
    x_ax = np.linspace(np.ceil(d2_min),np.floor(d2_max),int(np.floor(d2_max)-np.ceil(d2_min)+1.))
    x_ax=x_ax[0:-1:np.round(len(x_ax)/3)]
    y_ax = np.linspace(np.ceil(d1_min),np.floor(d1_max),int(np.floor(d1_max)-np.ceil(d1_min)+1.))
    y_ax=y_ax[0:-1:np.round(len(y_ax)/3)]
    ax.set_xticks(x_ax)
    ax.set_yticks(y_ax)
    ax.set_zticks([zmin,0,zmax])
    ax.xaxis.set_tick_params(pad=-4)
    ax.yaxis.set_tick_params(pad=-4)
    ax.zaxis.set_tick_params(pad=-2)
    
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    

    if zslice is not None:
        # Plane of DIP=0
        c_plane = colorConverter.to_rgba('k', alpha=0.2)
        verts = [array([(d1_min,d2_min), (d1_min,d2_max), (d1_max,d2_max), (d1_max,d2_min), (d1_min,d2_min)])]
        poly = PolyCollection(verts, facecolors=c_plane)
        ax.add_collection3d(poly, zs=[zslice], zdir='z')

        # Plot intersection of surface with DIP=0 plane
        CS = contour(np.log10(DD1),np.log10(DD2),E,levels=[zslice], linewidths=3, colors='k')

    # If needed, manually set zlim
    if zlim is not None: ax.set_zlim(zlim)
    
    if title is not None: ax.set_title(title)

    # Save plot, or show it
    if fname is None: show()
    else:
        plt.tight_layout()
        plt.savefig(fname)
        plt.clf()
        plt.close()


def plot_slices(d1_min, d1_max, d2_min, d2_max, E0, E1, E2, E3, h1, h2, r1,   \
                r1r, r2, r2r, alpha1, alpha2, scatter_points=None, N=50,      \
                fname=None, ylim=None, axes=None, drug1="drug1",              \
                drug2="drug2", units="M", zlabel="z", d1_slices=None,         \
                d2_slices=None, d1buffer=0., d2buffer=0., title=None):
    
    if axes is None or len(axes)!=2:
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
    else:
        ax1, ax2 = axes
    
    d1 = logspace(d1_min-d1buffer, d1_max, N)
    d2 = logspace(d2_min-d2buffer, d2_max, N)
    
    if d1_slices is None: d1_slices = [0, np.power(10.,d1_max)]
    else: d1_slices = [np.power(10.,i) for i in d1_slices]
    if d2_slices is None: d2_slices = [0, np.power(10.,d2_max)]
    else: d2_slices = [np.power(10.,i) for i in d2_slices]
    
    end_colors=["#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    
    for i,d in enumerate(d2_slices):
        E = hill_2D(d1, d, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
        if d < 0.01 and d > 0: label="%0.2e %s %s"%(d,units,drug2)
        else: label="%0.2f %s %s"%(d,units,drug2)
        lw=3
        if i == 0: c=end_colors[0]
        elif i == len(d2_slices)-1: c=end_colors[2]
        else:
            c="#555555"
            lw=1
            label="_nolabel_"
        
        ax1.plot(d1, E, c=c, linewidth=lw, label=label)
        if (i==0 or i==len(d2_slices)-1) and scatter_points is not None:
            sp = scatter_points.loc[scatter_points['drug2.conc']==d]
            ax1.scatter(sp['drug1.conc'], sp['effect'], c=c, edgecolors='k', s=20, label="_nolabel_", lw=0.5)
    
    for i,d in enumerate(d1_slices):
        E = hill_2D(d, d2, E0, E1, E2, E3, h1, h2, alpha1, alpha2, r1, r1r, r2, r2r)
        if d < 0.01 and d > 0: label="%0.2e %s %s"%(d,units,drug1)
        else: label="%0.2f %s %s"%(d,units,drug1)
        
        lw=3
        if i == 0: c=end_colors[1]
        elif i == len(d2_slices)-1: c=end_colors[3]
        else:
            c="#555555"
            lw=1
            label="_nolabel_"        
        ax2.plot(d2, E, c=c, linewidth=lw, label=label)
        if (i==0 or i==len(d1_slices)-1) and scatter_points is not None:
            sp = scatter_points.loc[scatter_points['drug1.conc']==d]
            ax2.scatter(sp['drug2.conc'], sp['effect'], c=c, edgecolors='k', s=20, label="_nolabel_", lw=0.5)
        
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    
    ax1.set_xlabel("%s (log[%s])"%(drug1,units))
    ax2.set_xlabel("%s (log[%s])"%(drug2,units))
    ax1.set_ylabel(zlabel)
    ax2.set_ylabel(zlabel)
    ax1.legend()
    ax2.legend()
    if ylim is not None:
        ax1.set_ylim(ylim)
        ax2.set_ylim(ylim)
    
    if axes is None:
        if title is not None: ax1.set_title(title)
        plt.tight_layout()
    
    if fname is None: show()
    else: plt.savefig(fname)
    clf()
    plt.close()

def plot_surface_contour(x, y, z, xlabel="x", ylabel="y", title="Plot", fname=None,
                         semilog=False, semilogx=False, semilogy=False, figsize=(4,3), 
                         cmap="inferno", center_cmap_on_zero=True, aspect='equal', 
                         vmin=None, vmax=None, crosshairs=False, power=1., levels=None):
    # Plot contours
    if (len(x.shape) == 1): XX, YY = np.meshgrid(x,y)
    else:
        XX, YY = x,y
        x = x[0,:]
        y = y[:,0]
    
    # These are tricks to make sure the colorbar in the legend comes out looking nice, even with this alpha
    alpha = 0.5    
    cmap = plt.get_cmap(cmap)(np.arange(256))
    
    # Edit the colormap to just show the redone colors
    for i in range(3): # Do not include the last column!
        cmap[:,i] = (1 - alpha) + alpha*cmap[:,i]
        cmap[:,i] = np.power(cmap[:,i],power)
    cmap = matplotlib.colors.ListedColormap(cmap, name='my_cmap')
    if (vmin is None and vmax is None):
        if center_cmap_on_zero:
            zmax = abs(np.nanmax(z))
            zmin = abs(np.nanmin(z))
            zmax = max(zmax,zmin)
            vmin = -zmax
            vmax = zmax
        else:
            vmin = np.nanmin(z)
            vmax = np.nanmax(z)
    cmap.set_bad(color="#AAAAAA")
    cmap.set_under('#AAAAAA')
                   
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, aspect=aspect)

    
    if levels is None:
        cs = contourf(XX,YY,z,cmap=cmap,vmin=vmin,vmax=vmax)
        CS = contour(XX, YY, z,cmap=cmap,vmin=vmin,vmax=vmax)
    else:
        cs = contourf(XX,YY,z,levels,cmap=cmap,vmin=vmin,vmax=vmax)
        CS = contour(XX, YY, z,levels,cmap=cmap,vmin=vmin,vmax=vmax)

    #ax.clabel(CS, inline=1, fontsize=8)
    ax.set_title(title)
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    # pcolor wants x and y to have dimensions 1 greater than z, so we have to resample x and y to match this dimensionality
    if semilog or semilogx:
        x_fudge = np.log10(x)
        xmin, xmax = min(x_fudge), max(x_fudge)
        x = np.logspace(xmin, xmax, x.shape[0]+1)
        ax.set_xscale('log')
    else:
        xmin, xmax = min(x), max(x)
        x = np.linspace(xmin, xmax, x.shape[0]+1)

    if semilog or semilogy:
        y_fudge = np.log10(y)
        ymin, ymax = min(y_fudge), max(y_fudge)
        y = np.logspace(ymin, ymax, y.shape[0]+1)
        ax.set_yscale('log')
    else:
        ymin, ymax = min(y), max(y)
        y = np.linspace(ymin, ymax, y.shape[0]+1)
    
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    
    ax.xaxis.set_tick_params(pad=-1)
    ax.yaxis.set_tick_params(pad=-1)


    
# 
#
    #pco = ax.pcolormesh(x, y, z, cmap=cmap, vmin=vmin, vmax=vmax,zorder=0)
#    
    #pco.set_edgecolor('face')

    #cbar = plt.colorbar(pco,ticks=[vmin, 0, vmax],orientation='horizontal')
    cbar = plt.colorbar(cs,orientation='horizontal')
    #cbar.ax.set_xticklabels(['%.2f'%vmin, 0, '%.2f'%vmax])
    
    cbar.set_ticks(cbar.get_ticks()[[0,int(floor(len(cbar.get_ticks())/2)),-1]])
    cbar.solids.set_rasterized(True)
    cbar.solids.set_edgecolor('face')
    cbar.outline.set_visible(False)

    
    if crosshairs:
        ax.plot([xmin,xmax],[0,0],'k--')
        ax.plot([1,1],[ymin,ymax],'k--')
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
    
    fig.tight_layout()
    if fname is None: plt.show()
    else:
        plt.savefig(fname)
        plt.cla()
        plt.clf()
        plt.close()
    return CS,ax
