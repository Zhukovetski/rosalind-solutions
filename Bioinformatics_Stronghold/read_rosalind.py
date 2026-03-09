import glob
import os


def read_latest_rosalind_file():
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

    with open(latest_file, mode="r", encoding="utf-8") as f:
        content = f.read().strip()

    return content
