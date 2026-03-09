from collections import Counter
import re

def solve_task(content):
    # Разделяем по заголовкам, игнорируем пустые элементы
    entries = [e.strip() for e in re.split(r"(?=>)", content) if e.strip()]
    
    max_gc = -1
    max_id = ""

    for entry in entries:
        # Разделяем на заголовок и тело
        lines = entry.splitlines()
        header = lines[0].replace(">", "")
        # Склеиваем ДНК, удаляя переносы строк
        sequence = "".join(lines[1:])
        
        counts = Counter(sequence)
        # Расчет процента (total теперь только по буквам A, T, G, C)
        gc_content = (counts['G'] + counts['C']) / len(sequence) * 100
        
        if gc_content > max_gc:
            max_gc = gc_content
            max_id = header

    return f"{max_id}\n{round(max_gc, 6)}"
