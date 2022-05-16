import math
import matplotlib.pyplot as plt

e = math.e
Y0 = -3  # Cauchy condition!


def f(x):
    return 3 * math.sqrt(x) / 2


def getPoints(N, a=0, b=1):
    h = (b - a) / N
    points = []
    cur_point = a
    for i in range(N + 1):
        points.append(cur_point)
        cur_point += h
    return points, h


def get_exact_solution(x):
    return x ** (3 / 2) - 3  # y = x^3/2 - 3


def get_euler_solution(points, h):
    result = [Y0]
    y = Y0

    for x in points[1:]:
        y = y + f(x) * h
        result.append(y)

    return result


def get_taylor_solution(points, h):
    f_dx1 = lambda x: 3 / (4 * math.sqrt(x))
    f_dx2 = lambda x: -3 / (8 * (x ** (3 / 2)))

    result = [Y0]
    y = Y0
    for x in points[1:]:
        y = y + f(x) * h + (h * h) / 2 * f_dx1(x) + (h * h * h) / 6 * f_dx2(x)
        result.append(y)

    return result


def runge_kutta(x, y, h):
    k1 = h * f(x)
    k2 = h * f(x + h / 2)
    k3 = h * f(x + h / 2)
    k4 = h * f(x + h)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6


def get_adams_solution(points, h):
    result = [Y0, runge_kutta(points[0], Y0, h)]
    y = Y0
    for x in points[1:-1]:
        y = y + h / 12 * (5 * f(x + h) + 8 * f(x) - f(x - h))
        result.append(y)
        x += h

    return result


def get_result(N):
    print()
    points, h = getPoints(N)
    print("N = " + str(N))
    print("h = ", h)
    print("values:", points)

    exact_solution = list(map(lambda x: get_exact_solution(x), points))
    print("Exact", exact_solution)

    euler_solution = get_euler_solution(points, h)
    print("Euler", euler_solution)

    taylor_solution = get_taylor_solution(points, h)
    print("Taylor K=3", taylor_solution)

    adams_solution = get_adams_solution(points, h)
    print("Adams K=2", adams_solution)

    plt.title(f'Exact and approximate solutions of F using {N} points')
    plt.plot(points, euler_solution)
    plt.plot(points, exact_solution)
    plt.plot(points, taylor_solution)
    plt.plot(points, adams_solution)
    plt.legend(["Exact", "Euler", "Taylor K=3", "Adams K=2"])
    plt.show()


if __name__ == '__main__':
    get_result(10)
    get_result(20)
    get_result(30)
