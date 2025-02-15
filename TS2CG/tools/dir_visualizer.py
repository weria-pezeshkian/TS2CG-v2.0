"""
CLI tool to place lipids to assign circular domains around inclusions or points.
"""

import argparse
from pathlib import Path
import numpy as np
import logging
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from ..core.point import Point

logger = logging.getLogger(__name__)


def make_fullscreen():
    """Make the Matplotlib figure full screen universally across platforms."""
    backend = matplotlib.get_backend()
    manager = plt.get_current_fig_manager()

    if backend == 'TkAgg':
        # For TkAgg backend (Linux, Windows, or macOS with Tkinter)
        manager.window.state('zoomed')  # Maximize window
    elif backend == 'Qt5Agg' or backend == 'QtAgg':
        # For Qt5Agg or QtAgg backend (requires PyQt5 or PyQt4)
        manager.window.showMaximized()
    elif backend == 'WXAgg':
        # For WXAgg backend
        manager.frame.Maximize(True)
    elif backend == 'MacOSX':
        # For MacOSX backend (native macOS backend)
        manager.full_screen_toggle()
    else:
        # Default: Manually set figure size to full screen
        fig = plt.gcf()
        fig.set_size_inches(16, 9)  # Adjust as needed for your resolution
        dpi = plt.rcParams['figure.dpi']
        fig.set_size_inches(1920 / dpi, 1080 / dpi)  # Example for 1920x1080

def draw_folder(membrane: Point, pointid: list, domain: bool, layer: str = "both",save=None,step=1,Proteins=False) -> None:
    """Assign lipids to domains based on curvature preferences"""
    layers = [membrane.outer]

    if membrane.monolayer or layer.lower() == "outer":
        layers = [membrane.outer]
        layer_names=["Outer"]
        layer_colors=["#0072B2"]
    elif layer.lower() == "inner":
        layers = [membrane.inner]
        layer_names=["Inner"]
        layer_colors=["#E69F00"]
    else:
        layers = [membrane.outer, membrane.inner]
        layer_names=["Outer","Inner"]
        layer_colors=["#0072B2","#E69F00"]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colorblind_1=["#009E73","#CC79A7","#56B4E9"]
    colorblind_2=['#332288', '#88CCEE', '#DDCC77', '#117733', '#882255','#CC6677', '#999999', '#F0E442', '#44AA99', '#AA4499','#D55E00', '#999933', '#66CCEE', '#882299', '#DD4499']
    second_run=False

    for i,layer in enumerate(layers):
        if not domain:
            x,y,z=layer.coordinates[:,0],layer.coordinates[:,1],layer.coordinates[:,2]
            ax.scatter(x[::step],y[::step],z[::step],s=50,label=layer_names[i],c=layer_colors[i])
        else:
            for k in range(0,16):
                reduced=layer.get_points_by_domain(domain_id=k)
                if not reduced.size == 0:
                    x,y,z=reduced[:,0],reduced[:,1],reduced[:,2]
                    if not second_run:
                        ax.scatter(x[::step],y[::step],z[::step],s=50,label=f"Domain {k}",c=colorblind_2[k])
                    else:
                        ax.scatter(x[::step],y[::step],z[::step],s=50,c=colorblind_2[k])


        labels=["c1","c2","c3"]
        for j,collection in enumerate(pointid):
            x=layer.coordinates[:,0][collection]
            y=layer.coordinates[:,1][collection]
            z=layer.coordinates[:,2][collection]
            if not x.size == 0:
                if not second_run:
                    ax.scatter(x,y,z,s=500,c=colorblind_1[j])
                else:
                    ax.scatter(x,y,z,s=500,label=labels[j],c=colorblind_1[j])
        second_run=True

    if Proteins:
        markers = ['D', 'd', 'p', '*', '+', 'x', 'h', 'H', '1', 'X']
        for i in range(len(markers)):
            ids=[Protein["point_id"] for Protein in membrane.inclusions.get_by_type(i)]
            if not len(ids) == 0:
                x=layers[0].coordinates[:,0][ids]
                y=layers[0].coordinates[:,1][ids]
                z=layers[0].coordinates[:,2][ids]
                ax.scatter(x,y,z,s=1000,label=f"Protein Type {i}",c="black",marker=markers[i])




    ax.legend(
    markerscale=1,  # Make legend markers larger
    labelspacing=1.5,  # Increase space between legend labels
    fontsize=16,  # Increase font size of legend labels
    handlelength=2,  # Make marker lines longer
    handleheight=2,  # Increase the height of the markers
    borderpad=1.5  # Increase padding around the legend
    )
    
    if save is not None:
        plt.savefig(save)
    else:
        make_fullscreen()
        plt.show()

def _get_centers(membrane,color1,color2,color3):
    colors=[color1,color2,color3]
    for i,item in enumerate(colors):
        colors[i]=list(map(int,list(filter(None,item.strip().split(",")))))
    return colors

def VIS(args) -> None:
    """Main entry point for Domain Placer tool"""
    parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p','--point_dir',default="point/",help="Specify the path to the point folder")
    parser.add_argument('-l','--leaflet',default="both",help="Choose which membrane leaflet to alter. Default is both")
    parser.add_argument('-S','--step',default=1,type=int,help="Only plot every nth point")
    parser.add_argument('-c1','--color1',default="",help="Highlights points in the visualization by point id. Excepts pointids like 3,7,22")
    parser.add_argument('-c2','--color2',default="",help="Highlights points in the visualization by point id. Excepts pointids like 3,7,22")
    parser.add_argument('-c3','--color3',default="",help="Highlights points in the visualization by point id. Excepts pointids like 3,7,22")
    parser.add_argument('-d','--Domain',default=False,action='store_true',help="Highlights points by domain, (up to 15 different domains are supported)")
    parser.add_argument('-P','--Protein',default=False,action='store_true',help="Shows proteins (up to 10 types are supported)")
    parser.add_argument('-s','--save_figure',default="",type=str,help="If path is outputfilename is provided, a figure will be saved")
    
   
    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO)

    if args.save_figure=="":
        filename=None
    else:
        filename=args.save_figure

    try:
        membrane = Point(args.point_dir)
        pointids=_get_centers(membrane,args.color1,args.color2,args.color3)
        draw_folder(membrane=membrane,layer=args.leaflet,pointid=pointids,domain=args.Domain,save=filename,step=args.step,Proteins=args.Protein)


    except Exception as e:
        logger.error(f"Error: {e}")
        raise
