import importlib.util
if importlib.util.find_spec("matplotlib") is not None:
    import matplotlib.pyplot as plt
import sys
import os


def read_input(path):
    """
    Чтение файла ввода.

    Args:
        path (str): путь до файла

    Returns:
        m (int): максимальное число точек смены высоты
        n (int): число точек, задающих кусочно-линейную функцию профиля рельефа
        x (list[float]): координаты рельефа по оси x
        y (list[float]): координаты рельефа по оси y
    """

    with open(path, 'r') as f:
        data = f.read().strip().split()

    if len(data) < 2:
        sys.exit("Недостаточно данных для построения высотной программы")
    if not data[0].isdigit() or int(data[0]) < 0:
        sys.exit("Некорректное число точек смены высоты")
    if not data[1].isdigit() or int(data[1]) < 0:
        sys.exit("Некорректное число точек рельефа")
    
    m, n = int(data[0]), int(data[1])

    if len(data[2:]) != 2 * n:
        sys.exit("Неправильный ввод - отсутсвует высота рельефа или количество точек не совпадает с N")

    x, y = [], []
    if n:
        x, y = list(map(float, data[2::2])), list(map(float, data[3::2]))
    return m, n, x, y


def write_output(path, change_points):
    """
    Запись в файл вывода.

    Args:
        path (str): путь до файла
        change_points (list[float]): список точек изменения высоты

    """

    with open(path, 'w') as f:
        f.write(str(len(change_points)) + '\n')
        if change_points:
            out = []
            for x, h in change_points:
                out.append(f"{x:.2f}")
                out.append(f"{h:.2f}")
            f.write(" ".join(out))


def add_points(x, y, n, n_seg = 50):
    """
    Функция добавления точек рельефа между точками изменения высоты.

    Args:
        x (list[float]): координаты рельефа по оси x
        y (list[float]): координаты рельефа по оси y
        n (int): число точек, задающих кусочно-линейную функцию профиля рельефа
        n_seg (int): количество сегментов, на которое делится отрезок рельефа, >=1

    Returns:
        x_added (list[float]): расширенный массив координат рельефа по оси x
        y_added (list[float]): расширенный массив координат рельефа по оси y
    """

    inclined_indexes = [i for i in range(len(y) - 1) if y[i] != y[i + 1]]  # массив индексов точек изменения высоты
    x_added, y_added = [], []

    for i in range(n):
        x_added.append(x[i])
        y_added.append(y[i]) 

        if i in inclined_indexes:
            dif_x = round((x[i + 1] - x[i]) / n_seg, 5)
            if dif_x < 0:
                sys.exit("Ошбка рельефа - дистанция должна увеличиваться")
            dif_y = round((y[i + 1] - y[i]) / n_seg, 5)
            for k in range(1, n_seg):    
                x_added.append(x[i] + dif_x * k)
                y_added.append(y[i] + dif_y * k)       

    return x_added, y_added


def math_err(x, y, i, j, h):
    """
    Функция подсчета площади(ошибки) между высотой и рельефом на отрезке [i, j] по оси x через 
    интеграл квадрата разности высоты полета и высоты рельефа.
        
    Args:
        x (list[float]): координаты рельефа по оси x
        y (list[float]): координаты рельефа по оси y 
        i (int): точка начального сегмента по оси x
        j (int): точка начального сегмента по оси x
        h (int): высота полета

    Returns:
        err (float): площадь(ошибка) сегмента
    """
    err = 0
    for k in range(i, j):
        err += round((x[k + 1] - x[k]) * (h**2 - h * y[k] - y[k + 1] * h + y[k + 1]*y[k] + ((y[k + 1] - y[k])**2) / 3), 5)
    return err


def merge_segments(x, y, m, plus_seg):
    """
    Функция построения оптимальной кусочно-постоянной траектории полёта над рельефом, 
    минимизирующей интеграл квадрата разности высоты полета и высоты рельефа (H - r(x))².

    Алгоритм:
        1. Сначала каждый участок [x[i], x[i+1]] рассматривается 
           как отдельный сегмент с высотой h = max(y[i], y[i+1]).
        2. Далее соседние сегменты постепенно объединяются попарно, 
           если это уменьшает или минимально увеличивает общую ошибку.
        3. Процесс продолжается, пока количество сегментов не станет равным 
           целевому значению m или меньше.                                        

    Args:
        x (list[float]): координаты рельефа по оси x
        y (list[float]): координаты рельефа по оси y
        m (int): максимальное число точек смены высоты
        plus_seg (int): число добавляемых фиктивных сегментов для учета возможости точек изменения высоты на краях рельефа

    Returns:
        segs (list[list[int, int, float, float]]): массив сегментов после объединения 
    """

    segs = []     
    for i in range(len(x) - 1):
        h = max(y[i], y[i + 1])
        err = math_err(x, y, i, i + 1, h)
        segs.append([i, i + 1, h, err])

    if y[0] == y[1]:
        segs[0][3] = -float("inf")
    # if y[-1] == y[-2]:
    #     segs[-1][3] = -float("inf")
    
    
    while len(segs) + plus_seg > m + 1:    
        best_k, best_dif_err, best_merge = None, float("inf"), None
        
        for k in range(len(segs) - 1):
            h_new = max(segs[k][2], segs[k + 1][2])
            err_new = math_err(x, y, segs[k][0], segs[k + 1][1], h_new)
            dif_err = err_new - (segs[k][3] + segs[k + 1][3]) 
            
            if dif_err < best_dif_err:
                best_k, best_dif_err = k, dif_err
                best_merge = [segs[k][0], segs[k + 1][1], h_new, err_new]
        
        if best_k != None:    
            segs[best_k] = best_merge
            segs.pop(best_k + 1)
        else:
            break

    return segs


def build_change_points(x, y, m, segs):
    """    
    Функция формирования из списка сегментов списка точек [x, h], в которых изменяется высота полёта БПЛА, 
    где x — координата по оси x(пройденное расстояние), h — высота полёта.
    Для начальной точки и каждой границы между соседними сегментами вычисляется высота полёта. 
    Если она отличается от предыдущей — добавляется точка изменения.                                    
           
    Args:
        x (list[float]): координаты рельефа по оси x
        y (list[float]): координаты рельефа по оси y
        m (int): максимальное число точек смены высоты
        segs (list[list[int, int, float, float]]: список сегментов после объединения

    Returns:
        change_points (list[tuple(float, float)]): список точек изменения высоты полета
    """

    change_points = []
    prev_h = y[0]

    for i in range(len(segs)):
        if abs(prev_h - segs[i][2]) > 1e-5:     
            change_points.append((x[segs[i][0]], segs[i][2]))
            prev_h = segs[i][2]

    if len(change_points) < m and abs(segs[-1][2] - y[-1]) > 1e-5:
        change_points.append((x[segs[-1][1]], y[segs[-1][1]]))

    return change_points


def plot_result(x, y, cp, x_added=None, y_added=None, title=""):
    plt.figure(figsize=(10, 5))
    plt.plot(x, y, "-o", label="Рельеф", linewidth=1)
    plt.plot(list(map(lambda x: x[0], cp)), list(map(lambda y: y[1], cp)), "ro", alpha=1, label="Точки изменения высоты")
    if x_added and y_added:
        plt.plot(x_added, y_added, ".", color="green", alpha=0.5, label="Добавленные точки")

    cp = [(x[0], y[0])] + cp 
    cp += [(x[-1], cp[-1][1])]
    points_x, points_y = [cp[0][0]], [cp[0][1]]
  
    for i in range(len(cp) - 1):
        points_x.extend([cp[i][0], cp[i + 1][0]])
        points_y.extend([cp[i][1], cp[i][1]])                  
    points_x.extend([cp[-2][0], cp[-1][0]])
    points_y.extend([cp[-2][1], cp[-1][1]])                

    plt.plot(points_x, points_y, label="Траектория БПЛА", linewidth=1)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlabel("Дальность, м")
    plt.ylabel("Высота, м")
    plt.title(title)
    plt.tight_layout()
    plt.show()


def compute_flight_plan(m, n, x, y, show_plot=False):
    """    
    Функция формирования из координат рельефа и максимального числа точек изменения высоты
    списка точек [x, h], в которых изменяется высота полёта БПЛА, 
    где x — координата по оси x (пройденное расстояние), h — высота полёта.                                    
           
    Args:
        m (int): максимальное число точек смены высоты
        n (int): количество координат рельефа
        x (list[float]): координаты рельефа по оси x
        y (list[float]): координаты рельефа по оси y
        show_plot (bool): флаг для отображения графика

    Returns:
        change_points (list[tuple(float, float)]): список точек изменения высоты полета
    """
     
    if (n <= 1) or (len(set(y)) == 1) or (m == 0 and y[0] == max(y)):  
        return []
    elif m == 0:
        sys.exit("Нельзя построить высотную программу по заданному рельефу.") 

    if m > 5000:
        num_add_points = 5000
    else:
        num_add_points = m

    x_added, y_added = add_points(x, y, n, num_add_points)

    plus_seg = int(y[0] < y[1])
    if n > 2 and y[-2] > y[-1]:
        plus_seg += 1
    # число добавляемых фиктивных сегментов при наличии наклонных краевых сегментах 
    # для учета возможости точек изменения высоты на краях рельефа

    segs = merge_segments(x_added, y_added, m, plus_seg)
    change_points = build_change_points(x_added, y_added, m, segs)
    if show_plot:
        try:
            plot_result(x, y, change_points, x_added, y_added, title="Высотная программа")
        except Exception as e:
            print(f"Не удалось построить график. Произошла ошибка {e}")
    return change_points 
    

def main():
    if len(sys.argv) < 3:
        sys.exit("Введите: python flight_program.py input.txt output.txt")
    inpath = sys.argv[1]
    outpath = sys.argv[2]
    do_plot = ("--plot" in sys.argv)

    if not os.path.exists(inpath):
        sys.exit("Некорректное имя входного файла")

    m, n, x, y = read_input(inpath)
    
    change_points = compute_flight_plan(m, n, x, y, do_plot)
    write_output(outpath, change_points)
    

if __name__ == "__main__":
    main()
    sys.exit(0)
