import DNA


class Indexer:
    def __init__(self, genome, k=11):
        self.index = {}
        self.k = k

        if isinstance(genome, DNA.DNA):
            for i in range(len(genome) - k + 1):
                kmer = str(genome[i : i + k])
                if kmer not in self.index:
                    self.index[kmer] = []
                if len(self.index.get(kmer)) < 1000:
                    self.index[kmer].append(i)

        else:
            kmer = ""
            for kmer, i in genome:
                if not kmer.ATGC:
                    continue
                kmer = str(kmer)
                if kmer not in self.index:
                    self.index[kmer] = []
                if len(self.index.get(kmer)) < 1000:
                    self.index[kmer].append(i)

            self.k = len(kmer)

    def query(self, kmer):
        return self.index.get(kmer, [])

    def add(self, kmer, pos):
        if len(self.index[kmer]) > 1000:
            print(f"{kmer} Слишком частый k-mer")
            return
        self.index[kmer].append(pos)


if __name__ == "__main__":
    # # 1. Создаем маленькую ДНК с известными повторами
    # # "ACGT" встречается на позиции 0 и 4
    # # "ZZZZ" - мусор
    # test_seq = DNA.DNA("ACGTACGTZZZZACGT")

    # # 2. Строим индекс с k=4
    # print("Строим индекс...")
    # indexer = Indexer(test_seq, k=4)

    # # 3. Проверяем то, что должно быть найдено
    # query = "ACGT"
    # result = indexer.query(query)
    # print(f"Поиск '{query}': {result}")

    # # Ожидание: [0, 4, 12]
    # # Реальность (с твоим кодом): Скорее всего [12] или ошибка.

    # # 4. Проверяем то, чего нет
    # print(f"Поиск 'AAAA': {indexer.query('AAAA')}")

    path = "/home/valentin/Downloads/ncbi_dataset/ncbi_dataset/data/GCA_000005845.2/GCA_000005845.2_ASM584v2_genomic.fna"
    # #path = "/home/valentin/chr1.fna"
    # name = "dna_fragment.csv"
    # substring = DNA.DNA("GGCTGGGCGTGGTGGCTCAC")
    # k = 11
    # distance = 500
    # dna = DNA.parse_fasta(path)
    # DNA.write_in_row(name, dna)
    # fragment = DNA.csv_read(name, fragment_size=k, offset=1)
    # indexer = Indexer(fragment)

    try:
        with open(path, "r") as f:
            lines = f.readlines()
            # Склеиваем всё, кроме заголовка
            raw_seq = "".join(
                line.strip() for line in lines if not line.startswith(">")
            )
            genome = DNA.DNA(raw_seq, "E.Coli")

        # Строим индекс
        idx = Indexer(genome, k=11)

        # ТЕСТ: Ищем кусок, который точно там есть (возьмем срез)
        real_pos = 10000
        query_seq = str(genome[real_pos : real_pos + 11])
        print(f"Looking for {query_seq} (Real pos: {real_pos})")

        matches = idx.query(query_seq)
        print(f"Found at: {matches}")

        if real_pos in matches:
            print("SUCCESS! Индекс работает.")
        else:
            print("FAILURE! Координата не найдена.")

    except FileNotFoundError:
        print("Файл не найден, проверь путь.")
