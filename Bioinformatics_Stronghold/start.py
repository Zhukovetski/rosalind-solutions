import glob
import os
import time
from itertools import tee

import matplotlib.pyplot as plt

import DNA


def get_latest_rosalind_file():
    """
    Ищет самый свежий файл rosalind_*.txt в папке загрузок.
    """
    # Получаем путь к папке Загрузки
    # Символ ~ автоматически заменяется на путь к домашней директории пользователя
    path_to_downloads = os.path.expanduser("~/Downloads")

    # 1. Формируем шаблон пути (например: C:/Users/Name/Downloads/rosalind_*.txt)
    # Звездочка * означает "любые символы"
    search_pattern = os.path.join(path_to_downloads, "rosalind_*.txt")

    # 2. Получаем список всех файлов, подходящих под шаблон
    list_of_files = glob.glob(search_pattern)

    # Если список пуст — возвращаем None
    if not list_of_files:
        return None

    # 3. Находим самый свежий файл
    # Функция max пройдется по списку файлов и сравнит их по критерию key
    latest_file = max(list_of_files, key=os.path.getctime)

    return latest_file


def plot_skew(dna_obj, fragment_size=10000):
    x = []
    y = []
    cumulative = 0

    # Мы будем считать "грубый" Cumulative Skew по окнам, это быстрее
    # Итерируемся по генератору, который ты уже написал!
    for pos, local_skew in dna_obj.skew():
        # local_skew - это (G-C)/Total.
        # Умножим на размер окна, чтобы получить чистую разницу (G-C) в буквах
        diff = local_skew * fragment_size

        cumulative += diff

        x.append(pos)
        y.append(cumulative)

    # Рисуем
    plt.figure(figsize=(15, 5))
    plt.plot(x, y, label="Cumulative Skew")
    plt.title(f"GC Skew Diagram for {dna_obj.name}")
    plt.xlabel("Position (bp)")
    plt.ylabel("Cumulative (G - C)")
    plt.grid(True)

    # Найдем минимум на графике
    min_y = min(y)
    min_x = x[y.index(min_y)]

    plt.axvline(min_x, color="r", linestyle="--", label=f"OriC approx ({min_x})")
    plt.legend()
    output_file = "gc_skew_plot.png"
    plt.savefig(output_file)
    print(f"График сохранен в файл: {output_file}")


# def solve_task(data):
#     lines = data.splitlines()
#     dna = DNA.DNAs(lines)
#     return dna[0].hamming_distance(dna[1])


# def main():


#     # Ищем файл
#     target_file = get_latest_rosalind_file()

#     if target_file:
#         print(f"Найден файл: {target_file}")

#         # Читаем файл
#         with open(target_file, "r") as f:
#             data = f.read().strip() # .strip() убирает лишние пробелы в начале и конце
#         print(data)
#         # Решаем задачу
#         answer = solve_task(data)

#         # Копируем ответ в буфер
#         pyperclip.copy(answer)
#         print(f"Ответ [{answer}] скопирован в буфер обмена!")

#         # ОПЦИОНАЛЬНО: Удаляем файл после использования, чтобы не мешался
#         # os.remove(target_file)
#         # print("Файл удален.")

#     else:
#         print("Файлы Rosalind не найдены в загрузках.")


def oric(path, fragment_size=10000, start=0, end=None):
    prew_skew = 0
    prew_position = start
    for dna in DNA.parse_fasta(path):
        for position, skew_value in dna.skew(fragment_size, start, end):
            if prew_skew * skew_value < 0:
                print(
                    f"OriC in {prew_position} - {min(position + fragment_size, len(dna))}"
                )
            prew_skew = skew_value
            prew_position = position

    # path = "/home/valentin/Downloads/ncbi_dataset/ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna"


if __name__ == "__main__":
    start = time.perf_counter()

    path = "/home/valentin/chr1.fna"
    name = "dna_fragment.csv"
    substring = DNA.DNA("GGCTGGGCGTGGTGGCTCAC")
    hm = 4
    distance = 500
    dna = DNA.parse_fasta(path)
    DNA.write_in_row(name, dna)
    fragment_size = 1000
    min_x = 0
    size = os.path.getsize(path)

    while fragment_size > 100:
        num_fragment = round(size**0.7)
        if num_fragment > 1000:
            num_fragment = 1000
        fragment = DNA.csv_read(
            name,
            min_x - fragment_size if min_x != 0 else 0,
            min_x + fragment_size if min_x != 0 else None,
            num_fragment,
        )
        gen1, gen2 = tee(fragment)
        first_element = next(gen1)
        fragment_size = first_element[2]
        min_x = DNA.pos_min_skew(gen2, fragment_size, image=True)
        size = fragment_size * 2

    print(str(min_x))
    # print(str(next(DNA.csv_read(name, min_x-50, min_x+50  + 9, 1))[1]))
    pos = DNA.DnaA(substring, name, min_x, hm, distance=distance)
    print(pos)
    pos = [p + min_x - distance for p in pos]
    print(pos)
    for i in range(len(pos)):
        print(str(next(DNA.csv_read(name, pos[i], pos[i] + len(substring), 1))[1]))

    finish = time.perf_counter()
    print(f"Время выполнения программы: {finish - start:.6f} секунд")


# for dna in DNA.parse_fasta(path):
#     plot_skew(dna, 250)


# oric(path)
# main()
# path = get_latest_rosalind_file()

# max_gc = 0
# for dna in DNA.parse_fasta("/home/valentin/Downloads/ncbi_dataset/ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna"):
#     gc = dna.gc_content()
# if dna.gc_content() > max_gc:
#     max_gc = gc
