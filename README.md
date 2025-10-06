# FASTA Parser

Python библиотека для работы с биологическими последовательностями в формате FASTA.

## Основные возможности

- **Класс Seq**: Работа с биологическими последовательностями (ДНК, РНК, белки)
- **Класс FastaReader**: Эффективное чтение и анализ FASTA файлов
- **Валидация данных**: Проверка формата и корректности последовательностей
- **Документация**: Полная HTML документация

## Установка

```bash
git clone https://github.com/username/fasta-parser.git
cd fasta-parser
pip install -r requirements.txt
```

## Быстрый старт

```python
from fasta_parser import Seq, FastaReader

# Работа с последовательностями
seq = Seq("ATGCGTAG", "Example DNA")
print(seq.get_alphabet_type())  # "DNA"
print(seq.get_gc_content())     # 37.5
print(len(seq))                 # 8

# Чтение FASTA файла
reader = FastaReader("sequences.fasta")
for sequence in reader:
    print(f">{sequence.header}")
    print(f"Type: {sequence.get_alphabet_type()}")
    print(f"Length: {len(sequence)}")
```

## Документация

Полная документация доступна в директории `docs/`. Для сборки HTML документации:

```bash
cd docs
pip install sphinx sphinx-rtd-theme
make html
```

Документация будет доступна в `docs/_build/html/index.html`

## Тестирование

Запуск всех тестов:

```bash
python -m unittest discover tests/
```

Запуск конкретного теста:

```bash
python -m unittest tests.test_seq
python -m unittest tests.test_fasta_reader
```

## Демонстрация

Запустите демонстрационную программу:

```bash
python examples/demo.py
```

## Структура проекта

```
fasta-parser/
├── fasta_parser/           # Основной пакет
│   ├── __init__.py
│   ├── seq.py             # Класс Seq
│   ├── fasta_reader.py    # Класс FastaReader
│   └── exceptions.py      # Пользовательские исключения
├── tests/                 # Модульные тесты
│   ├── test_seq.py
│   ├── test_fasta_reader.py
│   └── test_data/         # Тестовые FASTA файлы
│   └── demo.py
├── docs/                  # HTML документация
│   ├── conf.py
│   └── index.rst
├── README.md
├── requirements.txt
├── setup.py
└── LICENSE
```

## API Reference

### Класс Seq

Основные методы:
- `get_alphabet_type()` - определение типа последовательности
- `get_composition()` - анализ состава
- `reverse_complement()` - обратная комплементарность
- `translate()` - трансляция в белок

### Класс FastaReader

Основные методы:
- `read_sequences()` - генератор чтения последовательностей
- `get_file_stats()` - статистика файла
- `get_sequence_by_id()` - поиск по идентификатору
- `write_filtered_fasta()` - запись отфильтрованных данных


### Фильтрация последовательностей

```python
def filter_long_sequences(seq, min_length=100):
    return len(seq) >= min_length

reader = FastaReader("input.fasta")
count = reader.write_filtered_fasta(
    "long_sequences.fasta", 
    filter_func=filter_long_sequences,
    min_length=500
)
print(f"Сохранено {count} длинных последовательностей")
```
