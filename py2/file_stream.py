def read_2d_point(number) :
    li = []
    with open (f'.\..\{number}point.txt', 'r') as f :
        for line in f :
            xyPoint = line.split(',')
            li.append([])
            x = xyPoint[0]
            y = xyPoint[1][:-1]
            li[-1].append(float(x))
            li[-1].append(float(y))
    return li
    
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
    
def read_line(file_path) :
    """ read line from *.txt file.
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

def read_3d_point(number, state = True) :
    import math
    li = []
    if (state) :
        with open (f'.\..\{number}_3d_point.txt', 'r') as f :
            for line in f :
                xyPoint = line.split(',')
                li.append([])
                x = xyPoint[0]
                y = xyPoint[1]
                z = xyPoint[2][:-1]
                li[-1].append(float(x))
                li[-1].append(float(y))
                li[-1].append(float(z))
    else :
        with open (f'.\..\{number}_3d_point.txt', 'r') as f :
            for line in f :
                xyPoint = line.split(',')
                li.append([])
                
                latitude = float(xyPoint[0]) * math.pi / 180
                longitude = float(xyPoint[1]) * math.pi / 180
                altitude = float(xyPoint[2][:-1])
                 
                x = altitude * math.sin(latitude) * math.cos(longitude)
                y = altitude * math.sin(latitude) * math.sin(longitude)
                z = altitude * math.cos(latitude)
                
                li[-1].append(x)
                li[-1].append(y)
                li[-1].append(z)
        
    return li

def read_2d_answer(name) :
    li = []
    with open (f'.\..\{name}.txt', 'r') as f :
        for line in f :
            xyPoint = line.split(',')
            li.append([])
            x = xyPoint[0]
            y = xyPoint[1][:-1]
            li[-1].append(float(x))
            li[-1].append(float(y))
    return li

def read_3d_answer(name, state = True) :
    import math
    li = []
    if (state) :
        with open (f'.\..\{name}.txt', 'r') as f :
            for line in f :
                xyPoint = line.split(',')
                li.append([])
                x = xyPoint[0]
                y = xyPoint[1]
                z = xyPoint[2][:-1]
                li[-1].append(float(x))
                li[-1].append(float(y))
                li[-1].append(float(z))
    else :
        with open (f'.\..\{name}.txt', 'r') as f :
            for line in f :
                xyPoint = line.split(',')
                li.append([])
                
                latitude = float(xyPoint[0]) * math.pi / 180
                longitude = float(xyPoint[1]) * math.pi / 180
                altitude = float(xyPoint[2][:-1])
                 
                x = altitude * math.sin(latitude) * math.cos(longitude)
                y = altitude * math.sin(latitude) * math.sin(longitude)
                z = altitude * math.cos(latitude)
                
                li[-1].append(x)
                li[-1].append(y)
                li[-1].append(z)
    return li

def comparison_list(list1, list2) :
    for li in list1 :
        if li not in list2 :
            print('wrong')
            
            
if __name__ == '__main__' :
    std_answer_path = '3d_answer(all)'
    answer_path = '1000_3d_answer'
    
    std_answer = read_3d_answer(answer_path, False)
    answer = read_3d_answer(answer_path, False)