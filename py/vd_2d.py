import matplotlib.pyplot as plt

def read_point(file_path) :
    """ read point from *.txt file.

    Args:
        file_path (string): txt file path.

    Returns:
        list: return [[x,y]] coordinate list.
    """
    li = []
    with open (f'{file_path}', 'r') as f :
        first_run = True
        for line in f :
            if (first_run) :
                first_run = False
                continue
            xyPoint = line.split('\t')
            li.append([])
            x = xyPoint[1]
            y = xyPoint[2]
            li[-1].append(float(x))
            li[-1].append(float(y))
    return li

def read_vertex(file_path) :
    """ read vertex from *.txt file.

    Args:
        file_path (string): txt file path.

    Returns:
        list: return [[x,y]] coordinate list.
    """
    li = []
    with open (f'{file_path}', 'r') as f :
        for line in f :
            xyPoint = line.split(',')
            li.append([])
            x = xyPoint[0]
            y = xyPoint[1]
            li[-1].append(float(x))
            li[-1].append(float(y))
    return li

def read_edge(file_path) :
    """ read edge from *.txt file.

    Args:
        file_path (string): txt file path.

    Returns:
        list: return [[x,y]] coordinate list.
    """
    li = []
    with open (f'{file_path}', 'r') as f :
        for line in f :
            xyPoint = line.split(',')
            li.append([])
            start_x = xyPoint[0]
            start_y = xyPoint[1]
            end_x = xyPoint[2]
            end_y = xyPoint[3][:-1]
            li[-1].append(float(start_x))
            li[-1].append(float(start_y))
            li[-1].append(float(end_x))
            li[-1].append(float(end_y))
    return li

def read_face(name) :
    """ read face from *.txt file.

    Args:
        file_path (string): txt file path.

    Returns:
        list: return [[x,y]] coordinate list.
    """
    li = []
    with open (f'{name}', 'r') as f :
        for line in f :
            xyPoint = line.split(',')
            li.append([])
            x = xyPoint[0]
            y = xyPoint[1]
            li[-1].append(float(x))
            li[-1].append(float(y))
    return li

def draw_coordinate(_list, _c = 'black', _s = 5 , _alpha = 1, _marker = '.') :
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
    
def draw_arrow(_list, _color = 'black' ,_line_style = 'solid') :
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
        
def draw_line(_list, _color = 'black' ,_line_style = 'solid', _line_width = 0.5) :
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
        

import draw as draw
   
if __name__ == '__main__' :
    plt.axes().set_aspect('equal')  # xy 비율
    plt.xlim(-150000, 150000)
    plt.ylim(-150000, 150000)
    read_path = '.\\data\\point.txt'
    vertex_path = '.\\data\\answer_voronoi_vertex.txt'
    edge_path = '/home/jiminu/code/tsp-with-vd/data/answer_voronoi_edge.txt'
    face_path = '.\\data\\answer_voronoi_face.txt'
    
    temp_path = '/home/jiminu/code/tsp-with-vd/data/answer_voronoi_face.txt'
    
    # li = read_point(read_path)
    # vertices = read_vertex(vertex_path)
    edge = read_edge(edge_path)
    tem = read_edge(temp_path)
    # face = read_face(face_path)
    
    # draw_coordinate(li, _s = 30)
    # draw_coordinate(vertices, 'r', 5, _alpha = 0.5, _marker = 'v')
    draw_line(edge)
    draw_line(tem, 'green', _line_width = 2)
    # draw_line(unbounded_edge, 'red', 'dotted')
    # draw_arrow(edge)
    # draw_arrow(unbounded_edge, 'red', 'dotted')
    # draw_coordinate(face, 'b', 15, _alpha = 1)
    
    # queries_path = '.\\data\\answer_voronoi_queries.txt'
    # queries_edge = read_edge(queries_path)
    # draw_line(queries_edge, _line_width = 5)
    
    
    # ========== heuristic path ==========
    # result_file = 'tsp.txt'
    # info, result = draw.read_answer(result_file)
    # draw.draw_coordinate([face[int(result[0][0])]], 'red', 10)
    # draw.draw_coordinate([face[int(result[-1][0])]], 'red', 10)
    
    
    # result_line = []
    # for i in range(len(result)-1) :
    #     index = int(result[i][0])
    #     index2 = int(result[i+1][0])
    #     result_line.append([])
        # result_line[-1].append(face[index][0])
        # result_line[-1].append(face[index][1])
        # result_line[-1].append(face[index2][0])
        # result_line[-1].append(face[index2][1])
        
    # draw_line(result_line, 'green', _line_width = 2)
    # ====================================
    
    
    plt.savefig(f'data/vd.png')
    plt.show()