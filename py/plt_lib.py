from operator import length_hint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np

def make_3d_ax() :
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x_ratio = [[-10000], [10000]]
    y_ratio = [[-10000], [10000]]
    z_ratio = [[-10000], [10000]]
    ax.set_box_aspect((np.ptp(x_ratio), np.ptp(y_ratio), np.ptp(z_ratio)))
    
    return ax

def draw_2d_coordinate(_list, _c = 'black', _s = 5 , _alpha = 1, _marker = '.') :
    """ draw point from [[x,y]] coordinate list.

    Args:
        _list ([[x,y]]): [[x,y]] list.
        _c (str, optional): point color. Defaults to 'black'.
        _s (int, optional): point size. Defaults to 5.
        _alpha (int, optional): point transparency. Defaults to 1.
        _marker (str, optional): point style. Defaults to '.'.
    """
    xli = []
    yli = []
    for i in range(len(_list)) :
        xli.append(_list[i][0])
        yli.append(_list[i][1])
    plt.scatter(xli, yli, c= _c, s= _s, alpha= _alpha, marker= _marker)
    
def draw_3d_coordinate(ax, _list, _c = 'black', _s = 5, _alpha = 1, _marker = '.') :
    """ draw point from [[x,y,z]] coordinate list.

    Args:
        _list ([[x,y,z]]): [[x,y,z]] list.
        _c (str, optional): point color. Defaults to 'black'.
        _s (int, optional): point size. Defaults to 5.
        _alpha (int, optional): point transparency. Defaults to 1.
        _marker (str, optional): point style. Defaults to '.'.
    """
    xli = []
    yli = []
    zli = []
    for i in range(len(_list)) :
        xli.append(_list[i][0])
        yli.append(_list[i][1])
        zli.append(_list[i][2])
    ax.scatter(xli, yli, zli, c= _c, s= _s, alpha= _alpha, marker= _marker)
    
def draw_2d_arrow(_list, _color = 'black' ,_line_style = 'solid') :
    """ draw arrow.

    Args:
        _list ([[x,y]]): [[x,y]] list. x is start point, y is end point.
        _color (str, optional): arrow color. Defaults to 'black'.
        _line_style (str, optional): arrow style. Defaults to 'solid'.
    """
    for i in range(len(_list)) :
        start_x = _list[i][0]
        start_y = _list[i][1]
        end_x = _list[i][2]
        end_y = _list[i][3]
        plt.annotate('', xy=(end_x, end_y), xytext = (start_x, start_y),  arrowprops={'edgecolor': _color, 
                                                                                      'facecolor': 'red', 
                                                                                      'alpha' : 0.5, 
                                                                                      'arrowstyle' : '->',  
                                                                                      'linestyle': _line_style})
        
def draw_2d_line(_list, _color = 'black' ,_line_style = 'solid', _line_width = 0.5) :
    """ draw line.

    Args:
        _list ([[x,y]]): [[x,y]] list. x is start point, y is end point.
        _color (str, optional): line color. Defaults to 'black'.
        _line_style (str, optional): line style. Defaults to 'solid'.
        _line_width (float, optional): line width. Defaults to 0.5.
    """
    for i in range(len(_list)) :
        start_x = _list[i][0]
        start_y = _list[i][1]
        end_x = _list[i][2]
        end_y = _list[i][3]
        plt.plot([start_x, end_x], [start_y, end_y], color = _color, lw = _line_width)
        
    
def separate_xyz(li) :
    xli = []
    yli = []
    zli = []
    for i in li :
        xli.append(i[0])
        yli.append(i[1])
        zli.append(i[2])
    return xli, yli, zli
    
def dist_2d(p1, p2) :
    
    a1 = p1[0] - p2[0]
    b1 = p1[1] - p2[1]
    
    return math.sqrt(pow(a1,2) + (pow(b1,2)))

def dist_3d(p1, p2) :
    a1 = p1[0] - p2[0]
    b1 = p1[1] - p2[1]
    c1 = p1[2] - p2[2]

    return math.sqrt(pow(a1,2) + (pow(b1,2)) + (pow(c1, 2)))

def create_sphere(cx,cy,cz, r, resolution=10):
    phi = np.linspace(0, 2*np.pi, 2*resolution)
    theta = np.linspace(0, np.pi, resolution)

    theta, phi = np.meshgrid(theta, phi)

    r_xy = r*np.sin(theta)
    x = cx + np.cos(phi) * r_xy
    y = cy + np.sin(phi) * r_xy
    z = cz + r * np.cos(theta)

    return np.stack([x,y,z])
    
def draw_sphere_surface(ax, x, y, z, _color='black', _linewidth=0, _alpha=1, _rstride=1, _cstride=1) :
    ax.plot_surface(x, y, z, color=_color, linewidth=_linewidth, alpha=_alpha, rstride=_rstride, cstride=_cstride)
    
def draw_sphere_frame(ax, x, y, z, _color='black', _linewidth=0, _alpha=1, _rstride=1, _cstride=1) :
    ax.plot_wireframe(x, y, z, color=_color, linewidth=_linewidth, alpha=_alpha, rstride=_rstride, cstride=_cstride)

def draw_cuboid(ax, center, size, resolution = 8, _color = 'b', _linewidth = 1):
    # suppose axis direction: x: to left; y: to inside; z: to upper
    # get the (left, outside, bottom) point
    import numpy as np
    ox, oy, oz = center
    l, w, h = size

    x = np.linspace(ox-l/2,ox+l/2,num=resolution)
    y = np.linspace(oy-w/2,oy+w/2,num=resolution)
    z = np.linspace(oz-h/2,oz+h/2,num=resolution)
    x1, z1 = np.meshgrid(x, z)
    y11 = np.ones_like(x1)*(oy-w/2)
    y12 = np.ones_like(x1)*(oy+w/2)
    x2, y2 = np.meshgrid(x, y)
    z21 = np.ones_like(x2)*(oz-h/2)
    z22 = np.ones_like(x2)*(oz+h/2)
    y3, z3 = np.meshgrid(y, z)
    x31 = np.ones_like(y3)*(ox-l/2)
    x32 = np.ones_like(y3)*(ox+l/2)

    # outside surface
    ax.plot_wireframe(x1, y11, z1, color=_color, rstride=1, cstride=1, alpha=0.6, linewidth = _linewidth)
    # inside surface
    ax.plot_wireframe(x1, y12, z1, color=_color, rstride=1, cstride=1, alpha=0.6, linewidth = _linewidth)
    # bottom surface
    ax.plot_wireframe(x2, y2, z21, color=_color, rstride=1, cstride=1, alpha=0.6, linewidth = _linewidth)
    # upper surface
    ax.plot_wireframe(x2, y2, z22, color=_color, rstride=1, cstride=1, alpha=0.6, linewidth = _linewidth)
    # left surface
    ax.plot_wireframe(x31, y3, z3, color=_color, rstride=1, cstride=1, alpha=0.6, linewidth = _linewidth)
    # right surface
    ax.plot_wireframe(x32, y3, z3, color=_color, rstride=1, cstride=1, alpha=0.6, linewidth = _linewidth)

if __name__ == '__main__' :
    pass