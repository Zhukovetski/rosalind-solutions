import matplotlib.pyplot as plt
import numpy as np


def mandelbrot(h, w, max_iter=100):
    # Создаем сетку комплексных чисел c
    y, x = np.ogrid[-1.2 : 1.2 : h * 1j, -2 : 0.5 : w * 1j]
    c = x + y * 1j
    z = c
    # Массив для хранения итераций «побега»
    div_time = max_iter + np.zeros(z.shape, dtype=int)

    for i in range(max_iter):
        z = z**1.01 + c
        diverge = z * np.conj(z) > 4  # Проверка: |z| > 2
        escaping_now = diverge & (div_time == max_iter)
        div_time[escaping_now] = i  # Фиксируем итерацию побега
        z[diverge] = 2  # Предотвращаем переполнение (overflow)

    return div_time


plt.imshow(mandelbrot(800, 800), cmap="magma")
output_file = "mondelbrot1"
plt.savefig(output_file)
print(f"График сохранен в файл: {output_file}")
