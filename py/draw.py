import matplotlib.pyplot as plt
import numpy as np
import math

def read_point(name) :
    path = name
    li = []
    i = 0
    with open (path, 'r') as f :
        for line in f :
            if line == 'EOF' :
                break
            if i < 1 :
                i += 1
                continue
            xy = line.split('\t')
            li.append([])
            li[-1].append(float(xy[1]))
            li[-1].append(float(xy[2]))
    return li

def read_tsp_point(name) :
    path = name
    li = []
    i = 0
    with open (path, 'r') as f :
        for line in f :
            if line == 'EOF' :
                break
            if i < 1 :
                i += 1
                continue
            xy = line.split('\t')
            li.append([])
            li[-1].append(float(xy[1]))
            li[-1].append(float(xy[2]))
    return li

def read_answer(name) :
    path = name
    li = []
    info = []
    i = 0
    with open (path, 'r') as f :
        for line in f :
            if i < 8 :
                i += 1
                info.append(line)
                continue
            li.append([])
            x = line
            li[-1].append(int(x))
    return info, li

def draw_coordinate(_list, _c = 'black', _s = 5) :
    xli = []
    yli = []
    for i in range(len(_list)) :
        xli.append(_list[i][0])
        yli.append(_list[i][1])
    plt.scatter(xli, yli, c= _c, s= _s)
    
def dist(p1, p2) :
    
    a1 = p1[0] - p2[0]
    b1 = p1[1] - p2[1]
    
    return math.sqrt(pow(a1,2) + (pow(b1,2)))
    
if __name__ == '__main__' :
    num = 1000
    
    plt.axes().set_aspect('equal')  # xy 비율
    
    plt.grid(True, alpha=0.9)       # 그리드 투명도
    # plt.xticks(np.arange(5,100,10))   # 간격 설정 (시작, 끝, 간격)
    # plt.yticks(np.arange(5,100,10))
    # plt.gca().add_artist(plt.Circle([505,505], 105, fill = False))    # 원 그리기 (시작점, 길이, 채우기)
    
    li = read_point(num)
    answer = read_answer('answer')
    
    draw_coordinate(li)
    draw_coordinate([answer[0]], 'green', 20)
    plt.gca().add_artist(plt.Circle(answer[0], dist(answer[0], answer[-1]), fill = False, ec = 'green'))
    draw_coordinate(answer[1:], 'red', 10)
    plt.show()
    
    # for i in range(len(answer)) :
    #     if answer[i] != answer2[i] :
    #         print(f'{i} = {answer[i]} /// {answer2[i]}')