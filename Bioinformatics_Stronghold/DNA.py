import logging
import os
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

simple_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
detailed_formatter = logging.Formatter(
    "%(asctime)s [%(levelname)s] [Line: %(lineno)d] %(message)s"
)

# 2. Настраиваем обработчики
# Этот будет выводить ВСЁ в консоль (или основной файл)
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(simple_formatter)

# Этот будет записывать только WARNING и выше в отдельный файл с указанием строки КОДА
error_file_handler = logging.FileHandler("errors_debug.log")
error_file_handler.setLevel(logging.WARNING)
error_file_handler.setFormatter(detailed_formatter)

# 3. Регистрируем их
logging.basicConfig(
    level=logging.DEBUG,  # Ставим DEBUG, чтобы пропускать сообщения любого уровня к обработчикам
    handlers=[console_handler, error_file_handler],
)

logger = logging.getLogger(__name__)


class DNA:
    def __init__(self, sequence, name="Unknown"):

        self.sequence = sequence.upper()
        self.name = name
        self.ATCG = False

        valid_chars = set("ACGTURYSWKMBDHVN")
        current_chars = set(self.sequence)

        if not current_chars - set("ATCG"):
            self.ATCG = True

        invalid = set(self.sequence) - valid_chars
        if invalid:
            print(f"Warning: Strange chars {invalid} in {name}")

    def gc_content(self):
        # Расчет процента (total теперь только по буквам A, T, G, C)
        g_cnt = self.sequence.count("G")
        c_cnt = self.sequence.count("C")
        return (g_cnt + c_cnt) / len(self.sequence) * 100

    def counting_ACGT(self):
        # Количество нуклеотид ACGT в фрагмете ДНК
        count = Counter(self.sequence)
        return f"{count.get('A', 0)} {count.get('C', 0)} {count.get('G', 0)} {count.get('T', 0)}"

    def reverse_complement(self):
        # Вернуть реверсивную комплементарную ДНК
        trans_table = str.maketrans("ATCG", "TAGC")
        new_seq = self.sequence.translate(trans_table)[::-1]
        return DNA(new_seq)

    def transcribe(self):
        # Вернуть объект RNA (заменить T на U)
        return self.sequence.replace("T", "U")

    def hamming_distance(self, other):
        if not isinstance(other, DNA):
            raise TypeError("Can only add DNA to DNA")
        if len(self) != len(other):
            raise ValueError("Сan only be sequences of the same length")
        return sum(1 for c1, c2 in zip(self, other) if c1 != c2)

    def skew(self, end=None):
        if end is None:
            end = len(self)

        g_cnt = self.sequence.count("G")
        c_cnt = self.sequence.count("C")

        total = g_cnt + c_cnt
        skew_value = (g_cnt - c_cnt) / total if total > 0 else 0

        yield skew_value

    def substring(self, string, hm=0):
        pos = []
        lenght = len(string)
        for i in range(len(self.sequence) - lenght):
            sub = self[i : i + lenght]
            if sub.hamming_distance(string) <= hm:
                pos.append(i)

        return pos

    def __len__(self):
        return len(self.sequence)

    def __add__(self, other):
        if not isinstance(other, DNA):
            raise TypeError("Can only add DNA to DNA")
        return DNA(self.sequence + other.sequence)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return DNA(self.sequence[key])
        return self.sequence[key]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        # Красивый вывод: DNA(sequence='ATGC...')
        display_seq = self.sequence if len(self) < 15 else self.sequence[:12] + "..."
        return f"DNA(seq={display_seq}, len={len(self)})"


class DNAs:
    def __init__(self, sequence_list=None):
        self._storage = []
        if sequence_list is not None:
            for sequence in sequence_list:
                if isinstance(sequence, DNA):
                    self._storage.append(sequence)
                else:
                    self._storage.append(DNA(sequence))

    def __getitem__(self, index):
        return self._storage[index]

    def __repr__(self):
        return f"Vault({self._storage!r})"


def parse_fasta(filename, line_on=None):
    file_size = os.path.getsize(filename)
    current_id = None

    success_count = 0
    error_count = 0

    try:
        logger.info(f"Начало обработки файла: {filename}")
        with open(filename, "r") as f:
            with tqdm(
                total=file_size, unit="B", unit_scale=True, desc="Reading FASTA"
            ) as pbar:
                for line_num, line in enumerate(f, 1):
                    try:
                        pbar.update(len(line))

                        line = line.strip()
                        if not line:
                            continue

                        if line.startswith(">"):
                            # Сохраняем ПРЕДЫДУЩУЮ последовательность, если она есть
                            if current_id is not None:
                                success_count += 1
                                yield DNA(line, current_id)

                            # Начинаем собирать НОВУЮ
                            current_id = line[1:]
                        else:
                            if current_id is None:
                                logger.error(f"Строка {line_num}: Данные без заголовка")
                                error_count += 1
                                continue
                            yield DNA(line, current_id)

                    except ValueError as ve:
                        print(f"Ошибка валидации последовательности {current_id}: {ve}")
                        error_count += 1
                        continue

                if current_id is not None:
                    success_count += 1
                    yield DNA(line, current_id)

            logger.info(f"Завершено. Успешно: {success_count}, Ошибок: {error_count}")

    except FileNotFoundError:
        print(f"Критическая ошибка: Файл '{filename}' не найден.", file=sys.stderr)
    except PermissionError:
        print(
            f"Критическая ошибка: Нет прав доступа к файлу '{filename}'.",
            file=sys.stderr,
        )
    except Exception as e:
        print(f"Непредвиденная ошибка: {e}", file=sys.stderr)


def pos_min_skew(iterator, fragment_size=1, image=False):
    x = []
    y = []
    cumulative = 0
    for pos, dna, _ in iterator:
        for local_skew in dna.skew():
            cumulative += local_skew * fragment_size
            x.append(pos)
            y.append(cumulative)
    min_y = min(y)
    min_x = x[y.index(min_y)]
    if image:
        plt.figure(figsize=(15, 5))
        plt.plot(x, y, label="Cumulative Skew")
        plt.title("GC Skew Diagram")
        plt.xlabel("Position (bp)")
        plt.ylabel("Cumulative (G - C)")
        plt.grid(True)

        # Найдем минимум на графике
        min_y = min(y)
        min_x = x[y.index(min_y)]

        plt.axvline(min_x, color="r", linestyle="--", label=f"OriC approx ({min_x})")
        plt.legend()

        output_file = "gc_skew_plot.png"
        name, ext = os.path.splitext(output_file)
        counter = 1

        while os.path.exists(output_file):
            output_file = f"{name}({counter}){ext}"
            counter += 1
        plt.savefig(output_file)
        print(f"График сохранен в файл: {output_file}")
    return min_x


def csv_read(
    path, start=0, end=None, num_fragments=100, fragment_size=False, offset=False
):
    if end is None:
        end = os.path.getsize(path)
    size = end - start
    if not fragment_size:
        fragment_size = size // num_fragments
    elif not offset:
        num_fragments = size // fragment_size
    else:
        num_fragments = (size - fragment_size + 1) // offset
    with open(path, "r", encoding="ascii") as f:
        if not offset:
            for i in range(num_fragments):
                offset = start + i * fragment_size
                f.seek(offset)
                if i == num_fragments - 1:
                    chunk = f.read(end - offset)
                else:
                    chunk = f.read(fragment_size)
                yield offset, DNA(chunk), fragment_size

            else:
                for i in range(num_fragments):
                    f.seek(i)
                    chunk = f.read(fragment_size)
                    yield DNA(chunk), i


def write_in_row(name, dna):
    with open(name, "w", encoding="ascii") as f:
        for fragment in dna:
            f.write(str(fragment))


def DnaA(string, name, pos, hm=0, distance=500):
    dna = csv_read(name, pos - distance if pos > distance else 0, pos + distance, 1)
    fragment = next(dna)[1]
    pos = fragment.substring(string, hm)
    return pos


# region def Needleman_Wunsch(dna1, dna2):
#     len1, len2 = len(dna1), len(dna2)
#     dna1, dna2 = str(dna1), str(dna2)
#     score_matrix = np.zeros((len1 + 1, len2 + 1), dtype=int)

#     _match = 1
#     mismath = -1
#     gap = -1
#     score_matrix[0, :] = np.arange(len2 + 1) * gap
#     score_matrix[:, 0] = np.arange(len1 + 1) * gap

#     for i in range(1, len1 + 1):
#         for j in range(1, len2 + 1):
#             score_matrix[i, j] = max(
#                 score_matrix[i-1, j-1] +
#                 (_match if dna1[i-1] == dna2[j-1] else mismath),
#                 score_matrix[i-1, j] + gap,
#                 score_matrix[i, j-1] + gap)
#     print(score_matrix, end="\n\n")
#     result = np.full((3, max(len1, len2)+1), ' ', dtype="U1")

#     i, j = len1, len2
#     seq1, seq2, seqI = [], [], []
#     p_result = score_matrix.astype('U')
#     num = max(len1, len2) -1

#     while i > 0 or j > 0:
#         score = score_matrix[i, j]
#         score_diag = score_matrix[i-1, j-1]
#         if score == score_diag + (_match if dna1[i-1] == dna2[j-1] else mismath):
#             score = score_diag
#             p_result[i, j] = "\u2196"
#             i -= 1
#             j -= 1
#             result[0, num] = dna1[i]
#             result[1, num] = '|'
#             result[2, num] = dna2[j]
#             num -= 1
#             seq1.append(dna1[i])
#             seq2.append(dna2[j])
#             seqI.append("|")

#         elif score == score_matrix[i-1, j] + gap:
#             score = score_matrix[i-1, j]
#             p_result[i, j] = "\u2191"
#             i -= 1
#             result[0, num] = dna1[i]
#             result[1, num] = '|'
#             result[2, num] = '-'
#             num -= 1
#             seq1.append(dna1[i])
#             seq2.append('-')
#             seqI.append(' ')

#         else:
#             score = score_matrix[i, j-1]
#             p_result[i, j] = "\u2192"
#             j -= 1
#             result[0, num] = '-'
#             result[1, num] = '|'
#             result[2, num] = dna2[j]
#             num -= 1
#             seq2.append('-')
#             seq1.append(dna2[j])
#             seqI.append(' ')

#     print(pd.DataFrame(p_result).to_string(index=False, header=False), end="\n\n")
#     print(pd.DataFrame(result).to_string(index=False, header=False), end="\n\n")
#     seq1 = "".join(reversed(seq1))
#     seq2 = "".join(reversed(seq2))
#     seqI = "".join(reversed(seqI))
#     print(f"{seq1}\n{seqI}\n{seq2}")


# def Needleman_Wunsch(dna1, dna2):
#     dna1, dna2 = str(dna1), str(dna2)
#     i, j = -1, -1
#     i1, j1 = 0, 0
#     len1, len2 = len(dna1), len(dna2)
#     len_d = len2 - len1
#     i0 = j0 = gap = -3
#     match = 1
#     nums = {(-1, -1): 0}
#     one = (i, j)
#     two = ()
#     seq1, seq2 = [], []
#     num = max(i0 + gap, j0 + gap, match if dna1[i1] == dna2[j1] else -match)
#     nums[(0, 0)] = num

#     if nums[(i+1, j+1)] == nums[one] + match if dna1[i1] == dna2[j1] else -match:
#         i, j = i1, j1
#         seq1.append(dna1[i1])
#         seq2.append(dna2[j1])
#         one = (i, j)
#         i1 += 1
#         j1 += 1

#     elif j0 == nums[one] + gap:
#         j = j1
#         seq1.append(dna1[i1])
#         seq2.append("-")
#         one = (i, j)
#         ji += 1

#     elif i0 == nums[one] + gap:
#         i = i1
#         seq1.append("-")
#         seq2.append(dna2[j1])
#         one = (i, j)
#         i1 += 1


#     while i1 < len1 or j1 < len2:
#         if i1 - i:
#             n = i0
#             i0 = i0 + gap
#             for x in range(j1):
#                 num = max(i0, nums[(i, x)] + gap, match + n if dna1[i1] == dna2[x] else -match + n)
#                 n = nums[i, x]
#                 nums[(i1, x)] = num
#         if j1 - j:
#             m = j0
#             j0 = j0 + gap
#             for x in range(i1):
#                 num1 = max(nums[(x, j)]+ gap, j0, match + m if dna1[x] == dna2[j1] else -match + m)
#                 m = nums[x, j]
#                 nums[(x, j1)] = num1
#         s = nums[i1-1, j1-1]
#         num1 = max(nums[(i1, j)] + gap, nums[(i, j1)] + gap, match + s if dna1[i1] == dna2[j1] else -match + s)
#         nums[(i1, j1)] = num1
#         if nums[(i1, j1)] == nums[one] + match if dna1[i1] == dna2[j1] else -match:
#             if nums[(i, j1)] == nums[one] + gap:
#                 if i1 >= j1 - len_d:
#                     j = j1
#                     seq1.append("-")
#                     seq2.append(dna2[j1])
#                     one = (i, j)
#                     j1 += 1
#                 else:
#                     i = i1
#                     seq1.append(dna1[i1])
#                     seq2.append("-")
#                     one = (i, j)
#                     i1 += 1
#         else:
#             i, j = i1, j1
#             seq1.append(dna1[i1])
#             seq2.append(dna2[j1])
#             one = (i, j)
#             i1 += 1
#             j1 += 1


#     return "".join(seq1), "".join(seq2)
# endregion


def Needleman_Wunsch(dna1_obj, dna2_obj, Smith_Waterman=False):

    s1 = str(dna1_obj)
    s2 = str(dna2_obj)
    n, m = len(s1), len(s2)

    MATCH = 1
    MISMATCH = -1
    GAP = -1

    score_matrix = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(n + 1):
        score_matrix[i, 0] = i * GAP
    for j in range(m + 1):
        score_matrix[0, j] = j * GAP

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = score_matrix[i - 1, j - 1] + (
                MATCH if s1[i - 1] == s2[j - 1] else MISMATCH
            )
            delete_score = score_matrix[i - 1, j] + GAP
            insert_score = score_matrix[i, j - 1] + GAP

            score_matrix[i, j] = max(match_score, delete_score, insert_score)

    if Smith_Waterman:
        score_matrix[score_matrix < 0] = 0
        flat_index = np.argmax(score_matrix)
        i, j = np.unravel_index(flat_index, score_matrix.shape)
        n, m = i, j
    else:
        i, j = n, m

    align1, align2, symbol = [], [], []

    while i > 0 or j > 0:
        current_score = score_matrix[i, j]

        is_diag = False
        if i > 0 and j > 0:
            diagonal_score = score_matrix[i - 1, j - 1] + (
                MATCH if s1[i - 1] == s2[j - 1] else MISMATCH
            )
            if current_score == diagonal_score:
                is_diag = True

        is_up = False
        if i > 0:
            if current_score == score_matrix[i - 1, j] + GAP:
                is_up = True

        if is_diag:
            align1.append(s1[i - 1])
            align2.append(s2[j - 1])
            symbol.append("|")
            i -= 1
            j -= 1
        elif is_up:
            align1.append(s1[i - 1])
            align2.append("-")
            symbol.append(" ")
            i -= 1
        else:
            align1.append("-")
            align2.append(s2[j - 1])
            symbol.append(" ")
            j -= 1

        if Smith_Waterman and current_score == 0:
            break

    align1 = "".join(align1[::-1])
    align2 = "".join(align2[::-1])
    symbol = "".join(symbol[::-1])

    if Smith_Waterman:
        print("\nLocal Alignment Report:")
    else:
        print("\nGlobal Alignment Report:")
    print(f"Score: {score_matrix[n, m]}")
    print("-" * len(align1))
    print(align1)
    print(symbol)
    print(align2)
    print("-" * len(align1))


def Alignment_Algo(dna1_obj, dna2_obj, local=False):  # local=True -> Smith-Waterman
    s1 = str(dna1_obj)
    s2 = str(dna2_obj)
    n, m = len(s1), len(s2)

    MATCH = 1
    MISMATCH = -1
    GAP = -1

    score_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # 1. Инициализация
    # В локальном выравнивании края остаются нулями!
    if not local:
        for i in range(n + 1):
            score_matrix[i, 0] = i * GAP
        for j in range(m + 1):
            score_matrix[0, j] = j * GAP

    # 2. Заполнение
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = score_matrix[i - 1, j - 1] + (
                MATCH if s1[i - 1] == s2[j - 1] else MISMATCH
            )
            delete_score = score_matrix[i - 1, j] + GAP
            insert_score = score_matrix[i, j - 1] + GAP

            # Ключевой момент!
            if local:
                score_matrix[i, j] = max(0, match_score, delete_score, insert_score)
            else:
                score_matrix[i, j] = max(match_score, delete_score, insert_score)

    # 3. Выбор точки старта Traceback
    if local:
        # Начинаем с максимума во всей матрице
        flat_index = np.argmax(score_matrix)
        i, j = np.unravel_index(flat_index, score_matrix.shape)
    else:
        # Начинаем с конца
        i, j = n, m

    align1, align2, symbol = [], [], []

    # 4. Traceback
    # Добавили условие выхода для глобального (i>0 or j>0)
    # И условие выхода для локального (score > 0)
    while i > 0 or j > 0:
        current_score = score_matrix[i, j]

        # Для локального: если уперлись в 0 — выравнивание закончилось
        if local and current_score == 0:
            break

        # ВАЖНО: Приоритеты проверок. Сначала диагональ.
        is_diag = False
        if i > 0 and j > 0:
            diagonal_score = score_matrix[i - 1, j - 1] + (
                MATCH if s1[i - 1] == s2[j - 1] else MISMATCH
            )
            # В локальном мы могли прийти из 0, если диагональ дала плюс
            if current_score == diagonal_score or (
                local and current_score == 0 and diagonal_score < 0
            ):
                if current_score == diagonal_score:  # Строгая проверка для пути
                    is_diag = True

        is_up = False
        if i > 0:
            if current_score == score_matrix[i - 1, j] + GAP:
                is_up = True

        # Логика шага
        if is_diag:
            align1.append(s1[i - 1])
            align2.append(s2[j - 1])
            symbol.append("|" if s1[i - 1] == s2[j - 1] else " ")
            i -= 1
            j -= 1
        elif is_up:
            align1.append(s1[i - 1])
            align2.append("-")
            symbol.append(" ")
            i -= 1
        else:  # Left
            align1.append("-")
            align2.append(s2[j - 1])
            symbol.append(" ")
            j -= 1

    align1 = "".join(align1[::-1])
    align2 = "".join(align2[::-1])
    symbol = "".join(symbol[::-1])

    title = "Local (Smith-Waterman)" if local else "Global (Needleman-Wunsch)"
    print(f"\n{title} Report:")
    print(f"Max Score: {score_matrix.max() if local else score_matrix[n, m]}")
    print("-" * len(align1))
    print(align1)
    print(symbol)
    print(align2)
    print("-" * len(align1))


if __name__ == "__main__":
    d1 = DNA("ZZZZZATTACAZZZZZ")
    d2 = DNA("ATTACA")

    Alignment_Algo(d1, d2, local=True)

    # Needleman_Wunsch(d1, d2)
    # a, b = Needleman_Wunsch(dna1, dna2)
    # print(a)
    # print(b)
    # try:
    #     dna1 = DNA("ATGCGG")
    #     dna2 = DNA("AAAATTTT")

    #     print(f"DNA 1: {dna1}")
    #     print(f"GC Content: {dna1.gc_content():.2f}%")
    #     print(f"Counts: {dna1.counting_ACGT()}")

    #     # Магия сложения
    #     dna3 = dna1 + dna2
    #     print(f"Combined: {dna3}")

    #     # Reverse complement возвращает объект DNA
    #     rc = dna1.reverse_complement()
    #     print(f"RevComp: {rc}")
    #     print(dna1[:3])
    #     # Тест ошибки (раскомментируй, чтобы проверить)
    #     # bad_dna = DNA("ATGCZ")

    # except ValueError as e:
    #     print(f"Ошибка данных: {e}")
    # except Exception as e:
    #     print(f"Произошло что-то страшное: {e}")
